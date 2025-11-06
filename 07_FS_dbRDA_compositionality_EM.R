# -----------------------------------------------------------------------------#
# Forward selection and dbRDA of environmental variables and their impact 
# on EM communities, using the compositionally transformed community data 
# Original author: L. McKinley Nevins 
# Following: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
# April 24, 2025
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
#  aitchison distance matrix.                                                   #           
#                                                                               #
#################################################################################

####################
## 1. DATA PREP
##########@@########

# Load final phyloseq object of EM community that has transformed data  
ps_EM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final.RDS")

# Pull out OTU table for community analyses
spp <- otu_table(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

# Load in the aitchison distance matrix for the EM data 
load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_2025.Rdata")

# Load in the sample data file containing all environmental data 
env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_all.csv", row.names = 1)

# Trim down to just the columns with the environmental data in them 
env <- env[ -c(1:5) ]

#load original again for plotting later 
env2 <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_all.csv", row.names = 1)

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

#the threshold of R2 > 0.80 was used by Myers et al. 2013, I'm being a bit more conservative 

# Select a subset of environmental data that doesn't contain highly correlated variables 

env <- rownames_to_column(env, "Tree")

env_sub <- dplyr::select(env, Tree, lat, elev, mean_precip_mm, mean_summer_precip_mm, MAT, pct_N, 
                         K, B, ph, Sand, EC, avg_July_SPEI, count_mod_dry, count_sev_dry, apr1_SWE)

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

decorana(aitchison_EM) #length of DCA1 is 0.328, which I think means this is homogenous 


# Doing forward selection

# To use the ordiR2step function, you need to specify two models - the “empty” model 
#containing only intercept, and the “full” model including all variables (aka global 
#model). Also, the argument direction allows for selecting the selection direction 
#(forward, backward, both, with the default to both):

# Null model for stepwise selection
null_model <- vegan::capscale(aitchison_EM ~ 1, data = env_scaled)

# Full model
full_model <- vegan::capscale(aitchison_EM ~ ., data = env_scaled)

# Still getting a bit of redundancy in the drought and snow variables, but leaving them for now because I 
# am interested in these still as being biologically important 

# Forward selection
step_result <- ordistep(null_model, scope = formula(full_model), direction = "forward", permutations = 9999)

# View selected model
summary(step_result)


step_result
step_result$anova


# Sand, apr1_SWE, ph, lat, elev, count_sev_dry, MAT, and mean_precip_mm

# take env_scaled and merge host_ID to be able to model both 
hosts <- dplyr::select(env2, Host_ID)

mod_data <- merge(hosts, env_scaled, by = "row.names")

######################

#performing the dbRDA

#using just the 8 variables selected with forward selection plus host_ID
EM_RDA <- capscale(formula = aitchison_EM ~ Host_ID + Sand + apr1_SWE + ph + lat + elev + count_sev_dry + MAT + mean_precip_mm, 
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

# Only Host_ID significant 

summary(EM_RDA)

screeplot(EM_RDA)

EM_RDA

#                             Inertia            Proportion Rank
# Total         70.202745921404712703  1.000000000000000000     
# Constrained    3.853825820989005280  0.054895656436338619   15
# Unconstrained 66.348920100414460421  0.945104343563643590  404

# 5.5% of the variance is explained by the environmental variables and host 
# Inertia = total variation in the data 

#                              CAP1                CAP2                
# Eigenvalue            1.0610516464275297 0.52125777602509860 
# Proportion Explained  0.2753242351142981 0.13525722236489882
# Cumulative Proportion 0.2753242351142981 0.41058145747919694


#First axis explained 27.5% of the variation explained by the overall model (5.5%)


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
# enviornmental variables vs the host species 
env_df$Type <- c("Host", "Host", "Host", "Host", "Host", "Host", "Host", "Enviro", "Enviro", "Enviro",
          "Enviro", "Enviro", "Enviro", "Enviro", "Enviro")

env_df$Type <- as.factor(env_df$Type)

# get centroids for each site 
centroids <- site_df %>%
  group_by(Site) %>%
  summarize(CAP1 = mean(CAP1), CAP2 = mean(CAP2))


newSTorder <- c("Northern", "WFDP", "Andrews", "Southern")
site_df$Site <- factor(site_df$Site, levels = newSTorder)

# set colors for hosts 
                  # ABAM        ABGR        ABPR          ALRU         PSME         TABR         THPL       TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66")

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
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE")) +
  theme_bw() +
  labs(x = "CAP1 (27.53%)",
       y = "CAP2 (13.52%)",
       color = "Site") +
  coord_cartesian()  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

p

#######################################################

####################
## 2. DATA ANALYSIS
# PCA of forward selected variables
##########@@########


######## 

# Want to do some modeling of how the host identity and the environmental variables together contribute to 
# the community variation, and how much variation is explained by these factors 

## DATA PREP

### Run PCA to get simplified axes for these environmental variables of interest to use in modeling 

# Get subset of just the forward selected environmental variables 

env_PCA <- dplyr::select(env_sub, Sand, apr1_SWE, ph, lat, elev, count_sev_dry, MAT, mean_precip_mm) %>% rownames_to_column("Tree")

#PCA of Environmental Factors 
all.pca = prcomp(env_PCA[2:9], center = T, scale = T)

sd.all = all.pca$sdev
loadings.all = all.pca$rotation
env.names.all = colnames(env_PCA[2:9])
scores.all = as.data.frame(all.pca$x)
scores.all$Sample_ID = env2$Sample_ID
scores.all$Host_ID = env2$Host_ID
scores.all$Tree = env_PCA$Tree
summary(all.pca)


# PCA scores are 'scores.alltrees' with column for Sample_ID
loadings.all <- as.data.frame(loadings.all)

# get proportion of variance explained to add to each axis label 
pca_var <- all.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


# Plot the Results
PCA_plot <- ggplot(scores.all, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  geom_segment(data = loadings.all, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.all, aes(x = PC1 * 10, y = PC2 * 10, label = rownames(loadings.all)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_color_viridis_d(name = "Host_ID", option = "turbo") +
  labs(title = "PCA Biplot: Environmental Variation Across Host Tree Taxa",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host_ID") +
  theme(legend.position = "right")


PCA_plot

# Get top 4 factors for PC1
top_PC1 <- loadings.all[order(abs(loadings.all$PC1), decreasing = TRUE), ][1:3, ]

top_PC1 # count_sev_dry, elev, lat

# Get top 4 factors for PC2
top_PC2 <- loadings.all[order(abs(loadings.all$PC2), decreasing = TRUE), ][1:3, ]

top_PC2 # MAT, lat, elev

# Subset scores to get the position of each tree on the two PCA axes 

tree_scores <- dplyr::select(scores.all, PC1, PC2, Tree)


# Save this dataframe for us in the functional modeling analyses 

write.csv(tree_scores, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_tree_scores_PCs.csv")




##### SANDBOX ########




# Read in the distance to centroid values caclulated for the community compositionality analyses for Aitchison distance 
centroid_enviro_EM <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_centroid_distance_enviro.csv")

# Grab distance to centroid values 
centroid_enviro_EM <- dplyr::select(centroid_enviro_EM, Tree, Site, DistanceToCentroid)

# Grab scaled environmental variables of interest from above 
env_subset <- dplyr::select(env_scaled, Sand, apr1_SWE, ph, lat, elev, count_sev_dry, MAT, mean_precip_mm) %>% rownames_to_column("Tree")

# Load in sample data for the communities 
sample_metadata <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_all.csv")

sample_metadata$Tree <- sample_metadata$X

# add hosts that are missing from the distance to centroid data 
hosts <- dplyr::select(sample_metadata, Host_ID, Tree)


# Merge to make the dataset for modeling 
mod_data <- merge(centroid_enviro_EM, tree_scores, by = "Tree")
mod_data <- merge(mod_data, env_subset, by = "Tree")
mod_data <- merge(mod_data, hosts, by = "Tree")

hist(mod_data$Sand)
hist(mod_data$apr1_SWE)
hist(mod_data$ph)
hist(mod_data$lat)
hist(mod_data$elev)
hist(mod_data$count_sev_dry)
hist(mod_data$MAT)
hist(mod_data$mean_precip_mm)


## Dataset now has the scaled envionmental variables, as well as the position of each tree along the two PCA 
# axes of variation for those environmental variables 



# Run model with the environmental variables selected above
# Assessing the interaction of host identity with each of the environmental PCA axes 

start_time <- Sys.time()
mod1 <- lme4::lmer(DistanceToCentroid ~ Host_ID * (PC1 + PC2) + # Sand + apr1_SWE + ph + lat + elev + count_sev_dry + MAT + mean_precip_mm
                    ( 1| Site), # site as an additional random effect to account for other variation 
                   data    = mod_data, 
                   verbose = TRUE,
                   control = lmerControl(optimizer="bobyqa", mod.type = "lmer", optCtrl=list(maxfun=2e6)))
end_time <- Sys.time()


isSingular(mod1, tol = 1e-4) # FALSE

summary(mod1)
coef(mod1)
car::Anova(mod1, type=3)
performance::r2_nakagawa(mod1,
                         by_group   = FALSE,
                         tolerance  = 1e-05,
                         ci         = NULL,
                         iterations = 99)


## BOOTSTRAPPING ##

# Once each model has been run individually above, they can then be bootstrapped 

# This runs very quickly because the dataset is small 

#######
nboots <- 999

# Get fixed effect names
fixef_names <- names(fixef(mod1))

# Get random effect variance component names
randef_names <- names(VarCorr(mod1))
randef_names <- c(names(VarCorr(mod1)), "Residual")

# Create matrix to store both fixed and random effects
all_names <- c(fixef_names, randef_names)
boot_mat <- matrix(NA, nrow = nboots, ncol = length(all_names))
colnames(boot_mat) <- all_names

# ! ! ! ! TIMEWARN ! ! ! ! about XX minutes for 30 iterations ! ! ! ! ! ! ! ! ! ! ! ! ! !
for (b in 1:nboots) {
  cat('Bootstrap iteration', b, 'of', nboots, '....... \n')
  boot_data <- mod_data[sample(nrow(mod_data), replace = TRUE), ] # draw from dataset with replacement
  
  start_time <- Sys.time()
  tax_mod_boot <- try(lmer(DistanceToCentroid ~ Host_ID * (PC1 + PC2) +
                             ( 1| Site), # site as an additional random effect to account for other variation 
                      data    = boot_data, 
                      verbose = TRUE,
                      control = lmerControl(optimizer="bobyqa", mod.type = "lmer", 
                                            optCtrl=list(maxfun=2e9))))
  end_time <- Sys.time()
  
  
  # Skip if model failed
  if (inherits(tax_mod_boot, "try-error")) next
  
  # Fixed effects
  fixefs <- fixef(tax_mod_boot)
  
  # Random effects
  randefs <- as.data.frame(VarCorr(tax_mod_boot))
  randefs_vcov <- randefs$vcov
  names(randefs_vcov) <- randef_names
  
  # Combine
  all_effects <- c(fixefs, randefs_vcov)
  boot_mat[b, names(all_effects)] <- all_effects
  
  print(end_time - start_time)
}

#move the results array into a clearer name
tax_mod_boot999 <- as.data.frame(boot_mat)


isSingular(tax_mod_boot999, tol = 1e-4) # FALSE

summary(tax_mod_boot)
car::Anova(tax_mod_boot, type=3)
performance::r2_nakagawa(tax_mod_boot,
                         by_group   = FALSE,
                         tolerance  = 1e-05,
                         ci         = NULL,
                         iterations = 99)




#####Plotting Bootstrapped LMER Results #########################

# separate out the fixed and random effects from the bootstrapped results matrix, then get 
# the means and standard deviations of the effect sizes to plot 

#######

# right now done using the 30 bootstrapped results 

# Separate effects
fixed_df_tax <- dplyr::select(tax_mod_boot999, `(Intercept)`, PC1, PC2, Host_IDABGR, Host_IDABPR, Host_IDALRU, Host_IDPSME, Host_IDTABR, Host_IDTHPL, Host_IDTSHE,
                              `Host_IDABGR:PC1`, `Host_IDABPR:PC1`, `Host_IDALRU:PC1`, `Host_IDPSME:PC1`, `Host_IDTABR:PC1`, `Host_IDTHPL:PC1`, `Host_IDTSHE:PC1`,
                              `Host_IDABGR:PC2`, `Host_IDABPR:PC2`, `Host_IDALRU:PC2`, `Host_IDPSME:PC2`, `Host_IDTABR:PC2`, `Host_IDTHPL:PC2`, `Host_IDTSHE:PC2`) %>% as.data.frame() 


random_df_tax <- dplyr::select(tax_mod_boot999, Site, Residual) %>% as.data.frame()


# calculate means and standard errors 

# Initialize summary data.frame for fixed effects
fixed_summary_tax <- data.frame(
  Effect = colnames(fixed_df_tax),
  Mean   = apply(fixed_df_tax, 2, mean, na.rm = TRUE),
  SE     = apply(fixed_df_tax, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
)

# Compute confidence intervals from the SE
fixed_summary_tax$Lower_CI <- fixed_summary_tax$Mean - 1.96 * fixed_summary_tax$SE
fixed_summary_tax$Upper_CI <- fixed_summary_tax$Mean + 1.96 * fixed_summary_tax$SE

# Determine significance for each effect for plotting 
fixed_summary_tax$Significant <- with(fixed_summary_tax, Lower_CI > 0 | Upper_CI < 0)

# For random effects 
random_summary_tax <- data.frame(
  Effect = colnames(random_df_tax),
  Mean   = apply(random_df_tax, 2, mean, na.rm = TRUE),
  SE     = apply(random_df_tax, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
)

# Compute confidence intervals from the SE
random_summary_tax$Lower_CI <- random_summary_tax$Mean - 1.96 * random_summary_tax$SE
random_summary_tax$Upper_CI <- random_summary_tax$Mean + 1.96 * random_summary_tax$SE

# Determine significance for each effect for plotting 
random_summary_tax$Significant <- with(random_summary_tax, Lower_CI > 0 | Upper_CI < 0)


### PLOTTING EFFECT SIZES

# Get order of effect terms 
effect_names <- colnames(fixed_df_tax)  

# Convert to factor with levels in model order
fixed_summary_tax$Effect <- factor(fixed_summary_tax$Effect, levels = rev(effect_names))

fixed_plot_tax <- ggplot(fixed_summary_tax, aes(y = Effect, x = Mean, fill = Significant)) +
  geom_point(color='black', shape=21, size=3, aes(fill=factor(Significant)), show.legend = FALSE) + 
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.7) + 
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "red3", size=0.75) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  labs(y = "Fixed Effect", x = "Mean Effect Size ± 95% CI") +
  theme_minimal(base_size = 12)

fixed_plot_tax

## Random effects 
# Get order of effect terms 
rand_effect_names <- colnames(random_df_tax)  

# Convert to factor with levels in model order
random_summary_tax$Effect <- factor(random_summary_tax$Effect, levels = rev(rand_effect_names))

rand_plot_tax <- ggplot(random_summary_tax, aes(y = Effect, x = Mean, fill = Significant)) +
  geom_point(color='black', shape=21, size=3, aes(fill=factor(Significant)), show.legend = FALSE) + 
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.35) + 
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "red3", size=0.75) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  labs(y = "Random Effect", x = "Mean Effect Size ± 95% CI") +
  theme_minimal(base_size = 12)

rand_plot_tax


# pull together both plots to make a two-panel figure 
full_plot_tax <- fixed_plot_tax + rand_plot_tax + plot_layout(ncol = 2, widths = c(1, 1))

full_plot_tax


########################################################################################################





