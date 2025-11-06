# -----------------------------------------------------------------------------#
# Forward selection and dbRDA of environmental variables and their impact 
# on the functions of AM communities, using the compositionally transformed 
# community data 
# Original author: L. McKinley Nevins 
# Following: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
# April 24, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     phyloseq v 1.48.0
#                     vegan v 2.6.10
#                     adespatial v 0.3.24
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(adespatial); packageVersion("adespatial")
library(ggord): packageVersion("ggord")
library(ggordiplots): packageVersion("ggordiplots")
library(BiodiversityR); packageVersion("BiodiversityR")
library(ggvegan); packageVersion("ggvegan")
library(plotly); packageVersion("plotly")

#################################################################################
#                               Main workflow                                   #
#  Perform forward selection to narrow down the most important environmental    #
#  variables that explain variation in the community data. Then perform the     #
#  dbRDA ordination. Using the compositionally transformed data and the         #
#  aitchison distance matrix.                                                   #           
#                                                                               #
#################################################################################

# Load final phyloseq object of AM community that has transformed data  
ps_AM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")

# Pull out OTU table for community analyses
spp <- otu_table(ps_AM_clr) %>% as("matrix") %>% as.data.frame()

# Load in the aitchison distance matrix for the AM data 
load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_AM_traits.Rdata")

mat <- as.matrix(aitchison_AM_traits)

# 105 trees as many were excluded as not having functional variation different from the mean of the 
# community 


# Load in the sample data file containing the subset of environmental data 
env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_enviro_trait_aitchison_subset.csv", row.names = 1)

# Trim down to just the columns with the environmental data in them 
env <- env[ -c(1:5) ]

#load original again for plotting later 
env2 <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_enviro_trait_aitchison_subset.csv", row.names = 1)

#############################################################################

####evaluate the multicolinearity in the environmental data####
#ensure no variables are colinear before performing the forward selection 

cor_matrix <- cor(env, use = "everything")
print(cor_matrix)

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# Find variable pairs with absolute correlation > 0.70 (excluding self-correlation)
high_cor_pairs <- which(abs(cor_matrix) > 0.70 & upper.tri(cor_matrix), arr.ind = TRUE)

# Convert to a readable data frame
high_cor_df <- data.frame(
  Var1 = rownames(cor_matrix)[high_cor_pairs[,1]],
  Var2 = colnames(cor_matrix)[high_cor_pairs[,2]],
  Correlation = cor_matrix[high_cor_pairs]
)

print(high_cor_df)

#the threshold of R2 > 0.80 was used by Myers et al. 2013, I'm being a bit more conservative 

# Select a subset of environmental data that doesn't contain highly correlated variables 

env <- rownames_to_column(env, "Tree")

env_sub <- dplyr::select(env, Tree, elev, mean_precip_mm, MAT, pct_N, ph, Sand, count_mod_dry, apr1_SWE)

env_sub <- column_to_rownames(env_sub, "Tree")


# need to center and scale the environmental variables
env_scaled <- as.data.frame(scale(env_sub))


###########################################################################

### Forward Selection 

# Workflow follows closely: https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca

# No need to do the hellinger transform - since my matrix is already suitable for 
# euclidean methods 

# calculate DCA and consider the length of the first DCA axis to see how 
# homogenous or heterogenous the data are 
# This tells us if we can use a constrained method like RDA

decorana (aitchison_AM_traits) #length of DCA1 is 1.4383


# Doing forward selection

# To use the ordiR2step function, you need to specify two models - the “empty” model 
#containing only intercept, and the “full” model including all variables (aka global 
#model). Also, the argument direction allows for selecting the selection direction 
#(forward, backward, both, with the default to both):

# Null model for stepwise selection
null_model <- vegan::capscale(aitchison_AM_traits ~ 1, data = env_scaled)

# Full model
full_model <- vegan::capscale(aitchison_AM_traits ~ ., data = env_scaled)

# Forward selection
step_result <- ordistep(null_model, scope = formula(full_model), direction = "forward", permutations = 9999)

# View selected model
summary(step_result)


step_result
step_result$anova


# apr1_SWE, pct_N, count_mod_dry and elev significant 

# take env_scaled and merge host_ID to be able to model both 
hosts <- dplyr::select(env2, Host_ID)

mod_data <- merge(hosts, env_scaled, by = "row.names")

######################

#performing the dbRDA

#using just the 2 variables selected with forward selection 
AM_RDA <- capscale(formula = aitchison_AM_traits ~ Host_ID + apr1_SWE + pct_N + count_mod_dry + elev
                   , mod_data, distance = "euclidean", sqrt.dist = FALSE,
                   comm = NULL, add = FALSE, metaMDSdist = FALSE)


# p-value adjustment for multiple testing

# Run marginal tests for each environmental variable
anova_results <- anova(AM_RDA, by = "margin")

# Extract the raw p-values
raw_p <- anova_results$`Pr(>F)`

# Adjust p-value for multiple testing 
adjusted_p <- p.adjust(raw_p, method = "bonferroni")
# Add to the results table
anova_results$Adjusted_P <- adjusted_p
print(anova_results)

# Only pct_N and count_mod_dry still significant 

summary(AM_RDA)

screeplot(AM_RDA)

AM_RDA

# Partitioning of squared Euclidean distance:
#                 Inertia Proportion
# Total          132.24     1.0000
# Constrained     60.58     0.4581

# 45.8% of the variance is explained by the environmental variables 
# Inertia = total variation in the data 


# Importance of components:
#                          CAP1    CAP2      CAP3    
# Eigenvalue            59.0601 1.39298 0.1303568 
# Proportion Explained   0.4466 0.01053 0.0009857  
# Cumulative Proportion  0.4466 0.45713 0.4581187  

#First axis explained 45% of the variation explained by the overall model (45%)
#So the first axis is really carrying the team. 


plot(AM_RDA)

# Get site and biplot (environmental variable) scores
site_scores <- scores(AM_RDA, display = "sites")
env_scores  <- scores(AM_RDA, display = "bp")

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
env_df$Type <- c("Host", "Host", "Enviro", "Enviro", "Enviro", "Enviro")

env_df$Type <- as.factor(env_df$Type)

# get centroids for each site 
centroids <- site_df %>%
  group_by(Site) %>%
  summarize(CAP1 = mean(CAP1), CAP2 = mean(CAP2))


newSTorder <- c("Northern", "WFDP", "Andrews", "Southern")
site_df$Site <- factor(site_df$Site, levels = newSTorder)


#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# For just AM hosts 
AM_hosts <- c("#B93289FF", "#F48849FF", "#ffe24cFF")

sites <- c(15,16,17,18)


q <- ggplot() +
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
  scale_colour_manual(values=AM_hosts, 
                      name="Host Tree Species",
                      breaks=c("ALRU", "TABR", "THPL"),
                      labels=c("ALRU", "TABR", "THPL")) +
  theme_bw() +
  labs(x = "CAP1 (44.66%)",
       y = "CAP2 (10.53%)",
       color = "Site") +
  coord_cartesian()  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 11))

q



############

# INTERPRETATION

# This is explaining 45% of the variation and the environmental variables that are significant are 
# different from the taxonomic results. BUT, the shape of the ordination is nearly vertical, 
# which suggests that there are some issues with the data. 


