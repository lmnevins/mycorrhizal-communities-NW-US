# -----------------------------------------------------------------------------#
# Forward selection and dbRDA of environmental variables and their impact 
# on EM communities, using the compositionally transformed community data 
# Original author: L. McKinley Nevins 
# Following: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
# December 16, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     phyloseq v 1.48.0
#                     vegan v 2.6.10
#                     ggfortify v 0.4.17
#                     gginards v 0.2.0.1
#                     ggrepel v 0.9.6
#                     corrplot v 0.95
#                     car v 3.1.3
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(ggfortify); packageVersion("ggfortify")
require(gginnards); packageVersion("gginnards")
library(ggrepel); packageVersion("ggrepel")
library(corrplot); packageVersion("corrplot")
library(car); packageVersion("car")

#################################################################################
#                               Main workflow                                   #
#  Perform forward selection to narrow down the most important environmental    #
#  variables that explain variation in the community data. Then perform the     #
#  dbRDA ordination. Using the compositionally transformed data and the         #
#  aitchison distance matrix for the EM communities without THPL.               #           
#                                                                               #
#################################################################################

#################### -- 
## 1. DATA PREP
#################### -- 

# Load final phyloseq object of EM community that has transformed data  
ps_EM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final_no_THPL.RDS")

# Pull out OTU table for community analyses
spp <- otu_table(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

# Load in the aitchison distance matrix for the EM data 
load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_no_THPL.Rdata")

mat <- as.matrix(aitchison_EM)

# Load in the sample data file containing all environmental data 
env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_no_THPL.csv", row.names = 1)

# Trim down to just the columns with the environmental data in them 
env <- env[ -c(1:5) ]

#load original again for plotting later 
env2 <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_no_THPL.csv", row.names = 1)

#############################################################################

####evaluate the multicolinearity in the environmental data####
#ensure no variables are colinear before performing the forward selection 

cor_matrix <- cor(env, use = "everything")
print(cor_matrix)

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# Find variable pairs with absolute correlation > 0.7 (excluding self-correlation)
high_cor_pairs <- which(abs(cor_matrix) > 0.70 & upper.tri(cor_matrix), arr.ind = TRUE)

# Convert to a readable data frame
high_cor_df <- data.frame(
  Var1 = rownames(cor_matrix)[high_cor_pairs[,1]],
  Var2 = colnames(cor_matrix)[high_cor_pairs[,2]],
  Correlation = cor_matrix[high_cor_pairs]
)

print(high_cor_df)

#the threshold of R2 > 0.80 was used by Myers et al. 2013, I'm being a bit more conservative because the forward 
# selection process is more sensitive to auto-correlation.

# Select a subset of environmental data that doesn't contain highly correlated variables 

env <- rownames_to_column(env, "Tree")

env_sub <- dplyr::select(env, Tree, lat, elev, mean_precip_mm, mean_summer_precip_mm, MAT, pct_N, 
                         K, Fe, ph, Sand, EC, avg_July_SPEI, count_mod_dry)

env_sub <- column_to_rownames(env_sub, "Tree")


# need to center and scale the environmental variables
env_scaled <- as.data.frame(scale(env_sub))

###########################################################################

#####################
## 2. DATA ANALYSIS
# Forward Selection & dbRDA
####################

### Forward Selection 

# Workflow follows closely: https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca

# No need to do the hellinger transform - since my matrix is already suitable for 
# euclidean methods 

# calculate DCA and consider the length of the first DCA axis to see how 
# homogenous or heterogenous the data are 
# This tells us if we can use a constrained method like RDA

decorana(aitchison_EM) #length of DCA1 is 0.224, which I think means this is homogenous 


# Doing forward selection

# To use the ordiR2step function, you need to specify two models - the “empty” model 
#containing only intercept, and the “full” model including all variables (aka global 
#model). Also, the argument direction allows for selecting the selection direction 
#(forward, backward, both, with the default to both):

# Null model for stepwise selection
null_model <- vegan::capscale(aitchison_EM ~ 1, data = env_scaled)

# Full model
full_model <- vegan::capscale(aitchison_EM ~ ., data = env_scaled)

# Still getting a bit of redundancy in the drought variables, but leaving them for now because I 
# am interested in these still as being biologically important 

# Forward selection
step_result <- ordistep(null_model, scope = formula(full_model), direction = "forward", permutations = 9999)

# View selected model
summary(step_result)


step_result
step_result$anova


# pct_N + ph + lat + avg_July_SPEI + Sand + mean_precip_mm + count_mod_dry + K + EC + mean_summer_precip_mm + elev

# take env_scaled and merge host_ID to be able to model both 
hosts <- dplyr::select(env2, Host_ID)

mod_data <- merge(hosts, env_scaled, by = "row.names")

######################

#performing the dbRDA

#using just the 9 variables selected with forward selection plus host_ID
EM_RDA <- capscale(formula = aitchison_EM ~ Host_ID + pct_N + ph + lat + avg_July_SPEI + Sand + mean_precip_mm + 
                     count_mod_dry + K + EC + mean_summer_precip_mm + elev, 
                   mod_data, distance = "euclidean", sqrt.dist = FALSE,
                   comm = NULL, add = FALSE, metaMDSdist = FALSE)


# p-value adjustment for multiple testing

# Run marginal tests for each environmental variable
anova_results <- anova(EM_RDA, by = "margin")

# Extract the raw p-values
raw_p <- anova_results$`Pr(>F)`

# Adjust p-value for multiple testing 
adjusted_p <- p.adjust(raw_p, method = "bonferroni")
# Add to the results table
anova_results$Adjusted_P <- adjusted_p
print(anova_results)

# Only Host_ID (0.012), K (0.036), and EC (0.024) significant after p-value adjustment 

summary(EM_RDA)

screeplot(EM_RDA)

EM_RDA

# Partitioning of mean squared Euclidean distance:
#                  Inertia Proportion
# Total          3503.0    1.00000
# Constrained     309.3    0.08831
# Unconstrained  3193.6    0.91169

# 8.8% of the variance is explained by the environmental variables and host 
# Inertia = total variation in the data 

#                         CAP1    CAP2    
# Eigenvalue            71.7292 43.5389
# Proportion Explained   0.2319  0.1407 
# Cumulative Proportion  0.2319  0.3726


#First axis explained 23.2% of the variation explained by the overall model (8.8%)


plot(EM_RDA)

# Get site and biplot (environmental variable) scores
site_scores <- scores(EM_RDA, display = "sites")
env_scores  <- scores(EM_RDA, display = "bp")

# Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$Site <- env2$Site
site_df$Host_ID <- env2$Host_ID

env_df <- as.data.frame(env_scores)
env_df$Variable <- rownames(env_scores)

# edit the species names to be cleaner 
env_df$Variable <- gsub("Host_ID", "", env_df$Variable)


# Add a column that specifies the type of enviro data so I can change the font color for the 
# environmental variables vs the host species 
env_df$Type <- c("Host", "Host", "Host", "Host", "Host", "Host", "Enviro", "Enviro", "Enviro",
          "Enviro", "Enviro", "Enviro", "Enviro", "Enviro", "Enviro", "Enviro", "Enviro")

env_df$Type <- as.factor(env_df$Type)

# get centroids for each site 
centroids <- site_df %>%
  group_by(Site) %>%
  summarize(CAP1 = mean(CAP1), CAP2 = mean(CAP2))


newSTorder <- c("Northern", "WFDP", "Andrews", "Southern")
site_df$Site <- factor(site_df$Site, levels = newSTorder)

# set colors for hosts 
                  # ABAM        ABGR        ABPR          ALRU         PSME         TABR       TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#fffd66")

sites <- c(15,16,17,18)


p <- ggplot() +
  # Site points
  geom_point(data = site_df, aes(x = CAP1, y = CAP2, color = Host_ID, shape = Site), size = 2, alpha = 0.7) +
  # Arrows for environmental vectors
  geom_segment(data = env_df,
               aes(x = 0, y = 0, xend = CAP1 * 5, yend = CAP2 * 5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Text labels for environmental variables
  geom_text_repel(data = env_df,
                  aes(x = CAP1 * 5, y = CAP2 * 5, label = Variable, fontface = "bold"),
                  size = 4, color = "black") +
  scale_shape_manual(values=sites,
                     name="Site",
                     breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                     labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  scale_colour_manual(values=all_hosts, 
                      name="Host Tree Species",
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
  theme_bw() +
  labs(x = "CAP1 (23.19%)", # Remember to update to new values! 
       y = "CAP2 (14.07%)", # Remember to update to new values! 
       color = "Site") +
  coord_cartesian()  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

p

#######################################################

## -- END -- ## 
