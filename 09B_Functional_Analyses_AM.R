# -----------------------------------------------------------------------------#
# Chapter 2: Perform analyses on AM community functional composition
# Original Author: L. McKinley Nevins 
# December 16, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.2.1
#                     vegan v 2.6.10
#                     cluster v 2.1.8
#                     phyloseq v 1.48.0
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(vegan); packageVersion("vegan")
library(cluster); packageVersion("cluster")
library(phyloseq); packageVersion("phyloseq")

##################################################################################
#                               Main workflow                                    #
#  Load in AM functional variation summarized as the clr weighted representation #
#  of guilds across each site and host tree. Perform analyses of variation in    #
#  the functions to assess the drives of the variation.                          # 
#                                                                                #
##################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/"
setwd(wd)


# read in the clr transformed trait values per tree - for just the exploration types 

indiv_matrix <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_trait_clr_per_tree.csv")

# set sample_ID as rownames
indiv_matrix <- data.frame(indiv_matrix[,-1], row.names=indiv_matrix[,1])


# Read in metadata for each individual tree 

# Pull in phyloseq object with the clr transformed values 
ps_AM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")


indiv_meta <- data.frame(sample_data(ps_AM_clr))


#####################
## 3. DATA ANALYSIS
# PERMANOVA for traits x host 
# Beta dispersion of functions
####################

# calculate Aitchison distance using dist() from base R 

aitchison_AM_traits <- dist(indiv_matrix, method = "euclidean")

# Beta diversity of each tree in functional space 

# save Aitchison Distance matrix 
save(aitchison_AM_traits, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_AM_traits.Rdata")


## visualize ##

# PCA is appropriate for euclidean distances 
# Use the indiv_matrix with clr weighted trait values 

#### Perform PCA on the guilds for each individual tree 
traits_for_pca_AM <- merge(indiv_meta, indiv_matrix, by = 'row.names')


#PCA of differences in composition for Sites and Hosts
pca_traits_AM = prcomp(traits_for_pca_AM[36:39], center = T, scale = T)

sd.pca_traits_AM = pca_traits_AM$sdev
loadings.pca_traits_AM = pca_traits_AM$rotation
names.pca_traits_AM = colnames(traits_for_pca_AM[36:39])
scores.pca_traits_AM = as.data.frame(pca_traits_AM$x)
scores.pca_traits_AM$Site = traits_for_pca_AM$Site
scores.pca_traits_AM$Host_ID = traits_for_pca_AM$Host_ID
summary(pca_traits_AM)



# PCA scores are 'scores.pca_traits_AM' with column for Sites and Host_ID

loadings.pca_traits_AM <- as.data.frame(loadings.pca_traits_AM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_traits_AM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


scores.pca_traits_AM$Site <- as.factor(scores.pca_traits_AM$Site)
scores.pca_traits_AM$Host_ID <- as.factor(scores.pca_traits_AM$Host_ID)


#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# set colors for hosts 
#ALRU       TABR          THPL      
AM_hosts <- c("#B93289FF", "#F48849FF", "#ffe24cFF")

sites <- c(15,16,17,18)


# Plot the Results by site
PCA_traits_both_AM <- ggplot(scores.pca_traits_AM, aes(x = PC1, y = PC2, color = Host_ID, shape = Site)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_shape_manual(values=sites,
                     name="Site",
                     breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                     labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  scale_colour_manual(values=AM_hosts, 
                      name="Host Tree Species",
                      breaks=c("ALRU", "TABR", "THPL"),
                      labels=c("ALRU", "TABR", "THPL")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host_ID") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)) +
  guides(
    color = guide_legend(order = 1), # change legend order to be consistent across plots 
    shape = guide_legend(order = 2) 
  )

PCA_traits_both_AM


# Plot the Results by site alone
PCA_traits_site_AM <- ggplot(scores.pca_traits_AM, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Site), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Site") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)  
  )

PCA_traits_site_AM

# Plot the Results by host alone
PCA_traits_host_AM <- ggplot(scores.pca_traits_AM, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_colour_manual(values=AM_hosts, 
                      name="Host Tree Species",
                      breaks=c("ALRU", "TABR", "THPL"),
                      labels=c("ALRU", "TABR", "THPL")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Tree Species") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)  
  )

PCA_traits_host_AM


#####assess multivariate homogeneity of sites###########
#betadisper function in vegan 

#first object needs to be a dist object of Aitchison distances 
#second object is the groups of interest, as a vector
betadisper_traits_site <- vegan::betadisper(aitchison_AM_traits, indiv_meta$Site, type = "median", sqrt.dist = FALSE)

betadisper_traits_site


# Permutation test for homogeneity of multivariate dispersions
# p-value > 0.05, there is no evidence the groups have different dispersions
# p-value < 0.05, there is evidence the groups have different dispersions 

# This is permutational, so it will give a different p-value each time 
permutest(betadisper_traits_site)
#groups dispersions are not different

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      3   50.65  16.885 0.8992    999  0.449
# Residuals 126 2366.10  18.779



#if the dispersion is different between groups, then examine
scores(betadisper_traits_site, display = c("sites", "centroids"),
       choices = c(1,2))


#visualize 
plot(betadisper_traits_site, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C"))



boxplot(betadisper_traits_site)
mod.HSD <- TukeyHSD(betadisper_traits_site)
mod.HSD
plot(mod.HSD)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                          diff       lwr      upr     p adj
# Northern-Andrews  -0.4310349137 -3.033926 2.171857 0.9730192
# Southern-Andrews  -1.7442145022 -4.559721 1.071292 0.3750615
# WFDP-Andrews      -0.4318171487 -3.184516 2.320882 0.9768985
# Southern-Northern -1.3131795885 -4.252608 1.626249 0.6511106
# WFDP-Northern     -0.0007822351 -2.880108 2.878543 1.0000000
# WFDP-Southern      1.3123973535 -1.760475 4.385270 0.6828397


# Nicer boxplot 
distances_AM <- data.frame(
  Site = betadisper_traits_site$group,
  DistanceToCentroid = betadisper_traits_site$distances
)

# Define preferred order
site_order <- c("Northern", "WFDP", "Andrews", "Southern")

# Apply to your distances data
distances_AM$Site <- factor(distances_AM$Site, levels = site_order)

centroid_plot_AM <- ggplot(distances_AM, aes(x = Site, y = DistanceToCentroid, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal(base_size = 11) +
  labs(x = "Site",
       y = "Distance to Centroid") +
  scale_fill_manual(values=palette, 
                    name="Site",
                    breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                    labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

centroid_plot_AM


########


#check betadispersion between host tree taxa 
betadisper_traits_host <- betadisper(aitchison_AM_traits, indiv_meta$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper_traits_host

permutest(betadisper_traits_host)
#groups are different p = 0.001

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      2  227.37 113.686 8.9125    999  0.001 ***
#   Residuals 127 1620.00  12.756               


#if the dispersion is different between groups, then examine
plot(betadisper_traits_host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#B93289FF", "#F48849FF", "#ffe24cFF"))

boxplot(betadisper_traits_host)
mod.HSD_traits_host <- TukeyHSD(betadisper_traits_host)
mod.HSD_traits_host
plot(mod.HSD_traits_host)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#               diff        lwr        upr     p adj
# TABR-ALRU -3.290952 -5.1504784 -1.4314247 0.0001485
# THPL-ALRU -1.963897 -3.7739046 -0.1538889 0.0299956
# THPL-TABR  1.327055 -0.4714078  3.1255174 0.1908878

# TABR-ALRU = 0.0001 **
# THPL-ALRU = 0.03 *


# Nicer boxplot 
distances_AM2 <- data.frame(
  Host = betadisper_traits_host$group,
  DistanceToCentroid = betadisper_traits_host$distances
)

centroid_plot_AM2 <- ggplot(distances_AM2, aes(x = Host, y = DistanceToCentroid, fill = Host)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal(base_size = 11) +
  labs(x = "Host",
       y = "Distance to Centroid") +
  scale_fill_manual(values=AM_hosts, 
                    name="Host",
                    breaks=c("ALRU", "TABR", "THPL"),
                    labels=c("ALRU", "TABR", "THPL")) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

centroid_plot_AM2



# Run PERMANOVA models - for the individual tree level effects of both host and site 
# ** This is permutational so the results will be a little different every time. 

permanova.site <- vegan::adonis2(indiv_matrix ~ Site, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.site

# not significant, and only 3.8% explained


permanova.host <- adonis2(indiv_matrix ~ Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.host

# significant, and only 9.1% explained

permanova.both <- vegan::adonis2(indiv_matrix ~ Site * Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both

# significant, and 17.9% explained


permanova.both2 <- vegan::adonis2(indiv_matrix ~ Site + Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both2

# significant, and 13.5% explained


########################################################################################


