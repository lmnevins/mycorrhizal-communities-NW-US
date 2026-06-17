# -----------------------------------------------------------------------------#
# Chapter 2: Perform analyses on EM community functional composition
# Original Author: L. McKinley Nevins 
# December 16, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.2.1
#                     vegan v 2.6.10
#                     cluster v 2.1.8
#                     FD v 1.0.12.3
#                     phyloseq v 1.48.0
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(vegan); packageVersion("vegan")
library(cluster); packageVersion("cluster")
library(FD); packageVersion("FD")
library(phyloseq); packageVersion("phyloseq")
library(cowplot); packageVersion("cowplot")

##################################################################################
#                               Main workflow                                    #
#  Load in EM functional variation summarized as the clr weighted representation #
#  of exploration types across each site and host tree. Perform analyses of      #
#  variation in the functions to assess the drives of the variation.             # 
#                                                                                #
##################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/"
setwd(wd)


# read in the clr transformed trait values per tree - for just the exploration types 

indiv_matrix <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_trait_clr_per_tree.csv")

# set sample_ID as rownames
indiv_matrix <- data.frame(indiv_matrix[,-1], row.names=indiv_matrix[,1])



# Read in metadata for each individual tree 

# Pull in phyloseq object with the clr transformed values 
ps_EM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final_no_THPL.RDS")


indiv_meta <- data.frame(sample_data(ps_EM_clr))

#####################
## 2. DATA ANALYSIS
# PERMANOVA for traits x host 
# Beta dispersion of functions
##########@@########

# calculate Aitchison distance using dist() from base R 

aitchison_EM_traits <- dist(indiv_matrix, method = "euclidean")

mat <- aitchison_EM_traits %>% as.matrix() %>% as.data.frame()

# Beta diversity of each tree in functional space 

# save Aitchison Distance matrix 
save(aitchison_EM_traits, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_traits_no_THPL_2.0.Rdata")

## visualize ##

# PCA is appropriate for euclidean distances 
# Use the indiv_matrix with clr weighted trait values 

#### Perform PCA on the exploration type traits for each individual tree 
traits_for_pca_EM <- merge(indiv_meta, indiv_matrix, by = 'row.names')


#PCA of differences in composition for Sites and Hosts
# Scaling because these CLR values are not scaled 
pca_traits_EM = prcomp(traits_for_pca_EM[37:50], center = T, scale = T)

sd.pca_traits_EM = pca_traits_EM$sdev
loadings.pca_traits_EM = pca_traits_EM$rotation
names.pca_traits_EM = colnames(traits_for_pca_EM[37:50])
scores.pca_traits_EM = as.data.frame(pca_traits_EM$x)
scores.pca_traits_EM$Site = traits_for_pca_EM$Site
scores.pca_traits_EM$Host_ID = traits_for_pca_EM$Host_ID
summary(pca_traits_EM)



# PCA scores are 'scores.pca_traits_EM' with column for Sites and Host_ID

loadings.pca_traits_EM <- as.data.frame(loadings.pca_traits_EM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_traits_EM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


scores.pca_traits_EM$Site <- as.factor(scores.pca_traits_EM$Site)
scores.pca_traits_EM$Host_ID <- as.factor(scores.pca_traits_EM$Host_ID)


#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# set colors for hosts 
# ABAM        ABGR        ABPR          ALRU         PSME         TABR     TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#fffd66")

sites <- c(15,16,17,18)

# Shapes for all hosts 

# ABAM, ABGR, ABPR, ALRU, PSME, TABR, THPL, TSHE  
host_shapes <- c(15, 16, 17, 18, 7, 8, 9, 12)


# For just EM hosts 
EM_host_shapes <- c(15, 16, 17, 18, 7, 8, 12)


# Plot the Results by both site and host 
PCA_traits_both_EM <- ggplot(scores.pca_traits_EM, aes(x = PC1, y = PC2, color = Host_ID, shape = Site)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_shape_manual(values=sites,
                     name="Site",
                     breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                     labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  scale_colour_manual(values=all_hosts, 
                      name="Host Tree Species",
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host_ID") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

PCA_traits_both_EM


# Plot the Results by site alone
PCA_traits_site_EM <- ggplot(scores.pca_traits_EM, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Site), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 14) +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Site") +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 14),
        axis.text.y = element_text(colour="black", size = 14)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)) + 
  theme(legend.position = "none")

PCA_traits_site_EM


# Plot the Results by host alone
PCA_traits_host_EM <- ggplot(scores.pca_traits_EM, aes(x = PC1, y = PC2, color = Host_ID, shape = Host_ID)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 14) +
  scale_colour_manual(values=all_hosts, 
                      name="Host Tree Species",
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
  scale_shape_manual(values=EM_host_shapes, 
                     name="Host Tree Species",
                     breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                     labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Tree Species", shape = "Host Tree Species") +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 14),
        axis.text.y = element_text(colour="black", size = 14)) +
  theme(legend.position = "none")

PCA_traits_host_EM


#####assess multivariate homogeneity of sites###########
#betadisper function in vegan 

#first object needs to be a dist object of Aitchison distances 
#second object is the groups of interest, as a vector
betadisper_traits_site <- vegan::betadisper(aitchison_EM_traits, indiv_meta$Site, type = "median", sqrt.dist = FALSE)

betadisper_traits_site


# Permutation test for homogeneity of multivariate dispersions
# p-value > 0.05, there is no evidence the groups have different dispersions
# p-value < 0.05, there is evidence the groups have different dispersions 

# This is permutational, so it will give a different p-value each time 
permutest(betadisper_traits_site)
#groups dispersions are different 

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      3  448.3  149.44 12.894    999  0.001 ***
#   Residuals 367 4253.6   11.59                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


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
#                       diff        lwr        upr     p adj
# Northern-Andrews  -0.4265032 -1.6785164  0.8255100 0.8156702
# Southern-Andrews   1.0390434 -0.2443883  2.3224751 0.1584617
# WFDP-Andrews      -2.1120504 -3.4112871 -0.8128137 0.0001997 **
# Southern-Northern  1.4655466  0.1790187  2.7520744 0.0182571 *
# WFDP-Northern     -1.6855472 -2.9878425 -0.3832519 0.0050869 **
# WFDP-Southern     -3.1510938 -4.4836225 -1.8185650 0.0000000 ***

# Significant differences 

# WFDP-Andrews = 0.0002 **
# Southern-Northern = 0.018 *
# WFDP-Northern = 0.005 **
# WFDP-Southern = 0.0000000 ***


# Nicer boxplot 
distances_EM <- data.frame(
  Site = betadisper_traits_site$group,
  DistanceToCentroid = betadisper_traits_site$distances
)

# Define preferred order
site_order <- c("Northern", "WFDP", "Andrews", "Southern")

# Apply to distances data
distances_EM$Site <- factor(distances_EM$Site, levels = site_order)

centroid_plot_EM <- ggplot(distances_EM, aes(x = Site, y = DistanceToCentroid, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(x = "",
       y = "Distance to Centroid") +
  scale_fill_manual(values=palette, 
                    name="Site",
                    breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                    labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 14),
        axis.text.y = element_text(colour="black", size = 14)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)) + 
  theme(legend.position = "none")

centroid_plot_EM


########


# make dataframe of full species names 
sci_name <- c("A. amabilis", "A. grandis", "A. procera", "A. rubra", "P. menziesii", "T. brevifolia", "T. heterophylla")

Host_ID <- c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")


taxa <- data.frame(sci_name, Host_ID)

# Merge to all of the files

indiv_meta <- merge(indiv_meta, taxa, by = "Host_ID")


#check betadispersion between host tree taxa 
betadisper_traits_host <- betadisper(aitchison_EM_traits, indiv_meta$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper_traits_host

permutest(betadisper_traits_host)
#groups are different

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# # Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      6  359.6  59.935 5.1386    999  0.001 ***
# Residuals 364 4245.6  11.664                            
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 


#if the dispersion is different between groups, then examine
plot(betadisper_traits_host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#fffd66"))

boxplot(betadisper_traits_host)
mod.HSD_traits_host <- TukeyHSD(betadisper_traits_host)
mod.HSD_traits_host
plot(mod.HSD_traits_host)

#                diff        lwr         upr     p adj
# ABGR-ABAM  0.4464749 -1.4437815  2.33673128 0.9925316
# ABPR-ABAM  0.6608829 -1.3359491  2.65771494 0.9578171
# ALRU-ABAM -1.0918872 -3.1016438  0.91786946 0.6754566
# PSME-ABAM  1.1625870 -0.6939372  3.01911125 0.5106081
# TABR-ABAM -0.7482025 -2.6979051  1.20149996 0.9160225
# TSHE-ABAM -1.7504420 -3.6069662  0.10608220 0.0792599
# ABPR-ABGR  0.2144080 -1.8209516  2.24976770 0.9999244
# ALRU-ABGR -1.5383621 -3.5864032  0.50967909 0.2835810
# PSME-ABGR  0.7161122 -1.1817901  2.61401441 0.9221575
# TABR-ABGR -1.1946774 -3.1838206  0.79446578 0.5619727
# TSHE-ABGR -2.1969169 -4.0948191 -0.29901464 0.0117858 *
# ALRU-ABPR -1.7527701 -3.8995681  0.39402787 0.1928758
# PSME-ABPR  0.5017041 -1.5023672  2.50577546 0.9898082
# TABR-ABPR -1.4090855 -3.4997695  0.68159855 0.4174982
# TSHE-ABPR -2.4113249 -4.4153962 -0.40725360 0.0074101 **
# PSME-ALRU  2.2544742  0.2375247  4.27142378 0.0173760 *
# TABR-ALRU  0.3436846 -1.7593473  2.44671655 0.9990359
# TSHE-ALRU -0.6585548 -2.6755044  1.35839472 0.9604967
# TABR-PSME -1.9107896 -3.8679057  0.04632655 0.0607730
# TSHE-PSME -2.9130291 -4.7773375 -1.04872064 0.0001021 **
# TSHE-TABR -1.0022395 -2.9593556  0.95487668 0.7337081


# significant differences are 

# TSHE-ABGR = 0.011 
# TSHE-ABPR = 0.007
# PSME-ALRU = 0.017
# TSHE-PSME = 0.0001 


# Extract the p-values from the desired factor (column 4 of the output)
p_values <- mod.HSD_traits_host$group[, "p adj"]

# Load multcompView and get the compact letter display (CLD)
library(multcompView)
cld_letters <- multcompLetters(p_values)

# Extract just the letters
letters_df <- data.frame(Group = names(cld_letters$Letters), 
                         Letters = cld_letters$Letters)
print(letters_df)



# Nicer boxplot 
distances_EM2 <- data.frame(
  Host = betadisper_traits_host$group,
  DistanceToCentroid = betadisper_traits_host$distances
)


distances_EM2$Host_ID <- distances_EM2$Host

distances_EM2 <- merge(distances_EM2, taxa, by  = "Host_ID")


centroid_plot_EM2 <- ggplot(distances_EM2, aes(x = sci_name, y = DistanceToCentroid, fill = sci_name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(x = "",
       y = "Distance to Centroid") +
  scale_fill_manual(values=all_hosts, 
                    name="Host Tree Species",
                    breaks=c("A. amabilis", "A. grandis", "A. procera", "A. rubra", "P. menziesii", 
                             "T. brevifolia", "T. heterophylla"),
                    labels=c("A. amabilis", "A. grandis", "A. procera", "A. rubra", "P. menziesii", 
                             "T. brevifolia", "T. heterophylla")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 13, angle = 45, vjust = 0.5, face = "italic"),
        axis.text.y = element_text(colour="black", size = 14)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)) + 
  theme(legend.position = "none")

centroid_plot_EM2


# Run PERMANOVA models - for the individual tree level effects of both host and site 
# ** This is permutational so the results will be a little different every time. 

permanova.site <- vegan::adonis2(aitchison_EM_traits ~ Site, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.site

# significant, and only 4.9% explained


permanova.host <- adonis2(aitchison_EM_traits ~ Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.host

# significant, and only 5.0% explained

permanova.both <- vegan::adonis2(aitchison_EM_traits ~ Site * Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both

# significant, and only 19.0% explained


permanova.both2 <- vegan::adonis2(aitchison_EM_traits ~ Site + Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both2

# significant, and only 9.8% explained


########################################################################################

############################################## -- 
# (4) GATHER SIGNIFICANT RESULTS OF INTEREST  
############################################## -- 

# Gather up the two PCA plots and organize them for the figure

# Gather PCA plots:
# PCA_traits_host_EM, PCA_traits_site_EM

PCA_plots <- plot_grid(PCA_traits_site_EM, PCA_traits_host_EM, 
                       ncol = 1, nrow = 2, labels = c('(f)', '(h)'), label_size = 16, 
                       hjust = -0.2)

PCA_plots

ggsave("~/Dropbox/WSU/Mycorrhizae_Project/Publication_Materials/Figures/EM_trait_PCAs.png", 
       plot = PCA_plots, width = 5, height = 10, units = "in", dpi = 300)



# Gather up the two distance to centroid plots and organize them for the figure

# Gather distance to centroid plots:
# centroid_plot_EM, centroid_plot_EM2

centroid_plots <- plot_grid(centroid_plot_EM, centroid_plot_EM2,
                       ncol = 1, nrow = 2, labels = c('(f)', '(h)'), label_size = 16, 
                       hjust = -0.2)

centroid_plots

ggsave("~/Dropbox/WSU/Mycorrhizae_Project/Publication_Materials/Figures/EM_trait_centroids.png", 
       plot = centroid_plots, width = 5.25, height = 10, units = "in", dpi = 300)



## -- END -- ## 

