# -----------------------------------------------------------------------------#
# Transforming EM phyloseq data and doing compositional analyses
# Original Author: L. McKinley Nevins 
# December 14, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     phyloseq v 1.48.0
#                     purrr v 1.0.4
#                     Biostrings v 2.72.1
#                     corncob v 0.4.1
#                     vegan v 2.6.10
#                     patchwork v 1.3.1
#                     agricolae v 1.3.7
#                     rstatix v 0.7.2
#                     ggpubr v 0.6.0
#                     coin v 1.4.3
#                     compositions v 2.0.8
#                     multcompView v 0.1.10
#                     lme4 v 1.1.36
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(agricolae); packageVersion("agricolae")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(coin); packageVersion("coin")
library(compositions); packageVersion("compositions")
library(multcompView); packageVersion("multcompView")
library(lme4); packageVersion("lme4")

#################################################################################
#                               Main workflow                                   #
#  Transform data to account for compositionality, then perform analyses to     #
#  assess the taxonomic and functional composition of the communities, and to   #
#  compare community dissimilarity.                                             #
#                                                                               #
#################################################################################

# Load phyloseq object produced from the subset of the funguild community 
ps_EM <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_final_2025.RDS")

# Trim down just to include taxa that had functional assignments made 

# Load in taxon table for EM community
tax <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_subset_tax_2025.csv")

# pull out the genus and family levels, as well as the original ASV's. 
tax <- dplyr::select(tax, X, Family, Genus)

tax$Genus <- as.factor(tax$Genus)

summary(tax)


# Load in table of Genera with functional classifications 
exp <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/em_functions_no_THPL.csv")

exp$Genus <- as.factor(exp$Genus)

summary(exp)

# Merge datasets together 
matrix <- merge(tax, exp, by = "Genus", all.y = TRUE)

#get column name set for name merging down below
matrix$ASV <- matrix$X

# the number is correct - 1,627 ASVs with genus and functional assignments 

#start names file with the original ASVs 
ASV <- phyloseq::taxa_names(ps_EM)
ASV_long <- as.data.frame(ASV)

#change ASV names to something nicer to work with
taxa_names(ps_EM)
n_seqs <- seq(ntaxa(ps_EM))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_EM) <- paste("ASV", formatC(n_seqs, 
                                              width = len_n_seqs, 
                                              flag = "0"), sep = "_")
taxa_names(ps_EM)

# get shortened names 
ASV2 <- taxa_names(ps_EM) 
ASV_short <- as.data.frame(ASV2)

# join two dataframes
all_names <- cbind(ASV_long, ASV_short)

# save all_names file as key to the long and short ASV assignments 
write.csv(all_names, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/all_ASV_names_EM_2025.csv")


#merge taxa_names with the traits file to get the updated ASV names 
matrix <- merge(matrix, all_names, by = "ASV")


# Prune the final phyloseq object to just contain the taxa that have functional assignments 
taxa_to_keep <- matrix$ASV2

# Prune taxa from the phyloseq object that are NOT in the taxa_to_keep list 
ps_trimmed_funcs <- prune_taxa(taxa_to_keep, ps_EM)


# pull out values 
EM_final <- otu_table(ps_trimmed_funcs) %>% as("matrix") %>% as.data.frame()


# load this subset of taxa back into a phyloseq object 
ps_EM_final <- phyloseq(tax_table(tax_table(ps_EM)),
                      otu_table(otu_table(EM_final, taxa_are_rows=FALSE)),
                      sample_data(ps_EM))


#################################################################################################
# Getting to know your phyloseq data ####

# number of taxa - 1,627
ntaxa(ps_EM_final)

# number of samples - 372
nsamples(ps_EM_final)

# sample names
sample_names(ps_EM_final)
rank_names(ps_EM_final)
# taxa names
taxa_names(ps_EM_final)

# ASV table
otu_table(ps_EM_final) %>% View()

# how many sequences observed in each sample?
seq_counts <- otu_table(ps_EM_final) %>% rowSums() %>% as.data.frame()


####### 
# Note: There are 19 trees with <10 ASV's
######

# Remove one tree, S-ALRU-05, that doesn't have any reads 
ps_EM_final <- subset_samples(ps_EM_final, sample_names(ps_EM_final) != "S-ALRU-05_FWD")

# Check again
otu_table(ps_EM_final) %>% View() # Gone, now 371 trees 


# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_EM_final) %>% colSums()


# how many different samples was each taxon found in?
asv <- otu_table(ps_EM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

#inspect tax_table 

EM_tax_table <- as.data.frame(tax_table(ps_EM_final))

#count number of classifications in each column to determine coverage to taxonomic levels 

colSums(!is.na(EM_tax_table))

#  Kingdom  Phylum    Class    Order   Family   Genus    Species 
#  1627     1627      1627     1627    1627     1627     819
#  100%     100%      100%     100%    100%     100%     50.33%


# Load in sample data for the communities 
sample_metadata <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_no_THPL.csv")

sample_metadata <- column_to_rownames(sample_metadata, "X")


#############
# Save for untransformed final phyloseq 

saveRDS(ps_EM_final, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_func_final_no_THPL.RDS")

# This should be used for the differential abundance analyses later on 

#############

# Save 'matrix' which now has functional info matched up to OTU's and taxa
write.csv(matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_OTU_funcs_EM_no_THPL.csv")


####### more data exploration #########

#check rarefaction curve 
#using raw counts asv table 

otu.rare = otu_table(ps_EM_final)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

# we will use vegan rarecurve 
otu.rarecurve = rarecurve(otu.rare, step = 10000, label = F)

# very different sample depths 

#### 

# Pursuing a transformation because sequence data are inherently compositional, but standard rarefaction 
# is not repeatable and leads to the exclusion of a huge fraction of your data 

# Using centered log-ratio transformation per https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full

# Also: Van den Boogaart, K. G., and Tolosana-Delgado, R. (2013). Analyzing Compositional Data with R, 
# London, UK: Springer.

## Using robust clr transformation, which only transforms the non-zero values in the dataset, and is thus 
# less affected by high zero counts than standard clr is 

clr_EM <- decostand(asv, method = "rclr")

clr_EM <- as.data.frame(clr_EM)

# load the clr data back into a phyloseq object 
ps_EM_clr <- phyloseq(tax_table(tax_table(ps_EM_final)),
                     otu_table(otu_table(clr_EM, taxa_are_rows=FALSE)),
                     sample_data(sample_metadata))

#####################
##save phyloseq object with transformed data 
saveRDS(ps_EM_clr, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final_no_THPL.RDS")

## This is the phyloseq object to be used in all downstream analyses 
#####################

#####################Beta-diversity##########################################

## Beta-diversity can be calculated using Aitchison Distance, which is essentially the euclidean
# distance calculated between pairs of samples that have been transformed by CLR

# calculate Aitchison distance using dist() from base R 

aitchison_EM <- dist(clr_EM, method = "euclidean")

# save Aitchison Distance matrix 
save(aitchison_EM, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_no_THPL.Rdata")


## visualize ##

# PCA is appropriate for euclidean distances 

# get a few site-level columns to load back into the clr dataframe 

sample_data <-  data.frame(sample_data(ps_EM_clr))

sites <- dplyr::select(sample_data, Sample_ID, Site, Host_ID, Field_ID, Location)

clr_sites_EM <- merge(sites, clr_EM, by = 'row.names')

#PCA of differences in composition for Sites 
pca_clr_EM = prcomp(clr_sites_EM[7:1633], center = T, scale = F)

sd.pca_clr_EM = pca_clr_EM$sdev
loadings.pca_clr_EM = pca_clr_EM$rotation
names.pca_clr_EM = colnames(clr_sites_EM[7:1633])
scores.pca_clr_EM = as.data.frame(pca_clr_EM$x)
scores.pca_clr_EM$Site = clr_sites_EM$Site
scores.pca_clr_EM$Host_ID = clr_sites_EM$Host_ID
summary(pca_clr_EM)



# PCA scores are 'scores.pca_clr_EM' with column for Sites

loadings.pca_clr_EM <- as.data.frame(loadings.pca_clr_EM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_clr_EM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

#set colors for hosts
                 # ABAM        ABGR        ABPR          ALRU         PSME         TABR        TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#fffd66")

sites <- c(15,16,17,18)

# Plot the Results by site and host 
PCA_plot_both <- ggplot(scores.pca_clr_EM, aes(x = PC1, y = PC2, color = Host_ID, shape = Site )) +
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

PCA_plot_both 


# Plot the Results by site alone
PCA_plot_site <- ggplot(scores.pca_clr_EM, aes(x = PC1, y = PC2, color = Site)) +
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

PCA_plot_site

# Plot the Results by host alone
PCA_plot_host <- ggplot(scores.pca_clr_EM, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_colour_manual(values=all_hosts, 
                      name="Host Tree Species",
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
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

PCA_plot_host


#####assess multivariate homogeneity of sites###########
#betadisper function in vegan 

#first object needs to be a dist object of Aitchison distances 
#second object is the groups of interest, as a vector
betadisper.site <- vegan::betadisper(aitchison_EM, sample_data$Site, type = "median", sqrt.dist = FALSE)

betadisper.site

# Permutation test for homogeneity of multivariate dispersions
# p-value > 0.05, there is no evidence the groups have different dispersions
# p-value < 0.05, there is evidence the groups have different dispersions 

# This is permutational, so it will give a different p-value each time 
permutest(betadisper.site, pairwise = TRUE)
#groups dispersions are not significantly different (homogenous) p = 0.067
# Meets assumption for the permanova

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      3   69.6 23.2118 2.4404    999  0.068 .
# Residuals 367 3490.7  9.5114                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#if the dispersion is different between groups, then examine
scores(betadisper.site, display = c("sites", "centroids"),
       choices = c(1,2))


#visualize 
plot(betadisper.site, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C"))



boxplot(betadisper.site)
mod.HSD <- TukeyHSD(betadisper.site)
mod.HSD
plot(mod.HSD)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                        diff        lwr        upr     p adj
# Northern-Andrews   0.87035851 -0.2638289 2.00454596 0.1972249
# Southern-Andrews  -0.26121328 -1.4238624 0.90143588 0.9381006
# WFDP-Andrews       0.01073127 -1.1662355 1.18769807 0.9999953
# Southern-Northern -1.13157179 -2.2970258 0.03388217 0.0606610
# WFDP-Northern     -0.85962723 -2.0393648 0.32011032 0.2381911
# WFDP-Southern      0.27194456 -0.9351812 1.47907032 0.9376355

# No sig diffs between sites 

# Nicer boxplot 
distances_EM <- data.frame(
  Site = betadisper.site$group,
  DistanceToCentroid = betadisper.site$distances
)

# Define preferred order
site_order <- c("Northern", "WFDP", "Andrews", "Southern")

# Apply to your distances data
distances_EM$Site <- factor(distances_EM$Site, levels = site_order)

centroid_plot_EM <- ggplot(distances_EM, aes(x = Site, y = DistanceToCentroid, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal(base_size = 11) +
  labs(title = "",
       x = "",
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

centroid_plot_EM


permanova.site <- vegan::adonis2(aitchison_EM ~ Site, data = sample_data, method = "euclidean", permutations = 999)
permanova.site

# vegan::adonis2(formula = aitchison_EM ~ Site, data = sample_data, permutations = 999, method = "euclidean")
#              Df SumOfSqs      R2      F Pr(>F)    
# Model      3    333.1 0.01339 1.6599  0.001 ***
#   Residual 367  24546.9 0.98661                  
# Total    370  24879.9 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# adonis result is significant, and there was no heterogeneity of variances for the sites, so 
# this can be reported without concern 


permanova.both <- vegan::adonis2(aitchison_EM ~ Site * Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.both

# explained 9.6% 


permanova.both2 <- vegan::adonis2(aitchison_EM ~ Site + Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.both2

# explained 4%

# When including the interaction of tree host and site, 9.6% of the variation was explained. Environmental and host identity explain only 
# 9.6% of the variation. Only enviro explained 1%, ... just compare 

##################

#check betadispersion between host tree taxa 
betadisper.host <- betadisper(aitchison_EM, sample_data$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper.host

permutest(betadisper.host)
#groups are different p = 0.001
# Does not meet assumption for homogeneity of variances 

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#                Df  Sum Sq Mean Sq     F N.Perm Pr(>F)    
#   Groups      6  367.76  61.293 7.274    999  0.001 ***
#   Residuals 364 3067.19   8.426  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#if the dispersion is different between groups, then examine
plot(betadisper.host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#fffd66"))

boxplot(betadisper.host)
mod.HSD.host <- TukeyHSD(betadisper.host)
mod.HSD.host
plot(mod.HSD.host)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
#               diff        lwr         upr     p adj
# ABGR-ABAM  0.2630496 -1.3436091  1.86970825 0.9990257
# ABPR-ABAM -0.2231823 -1.9204269  1.47406231 0.9997240
# ALRU-ABAM -1.7368891 -3.4451192 -0.02865899 0.0433859
# PSME-ABAM  1.0349461 -0.5430413  2.61293345 0.4521910
# TABR-ABAM  1.5113030 -0.1458829  3.16848904 0.1001746
# TSHE-ABAM  1.3611215 -0.2168658  2.93910890 0.1424536
# ABPR-ABGR -0.4862319 -2.2162238  1.24376001 0.9813392
# ALRU-ABGR -1.9999387 -3.7407095 -0.25916793 0.0128427
# PSME-ABGR  0.7718965 -0.8412609  2.38505390 0.7915276
# TABR-ABGR  1.2482535 -0.4424559  2.93896281 0.3040096
# TSHE-ABGR  1.0980720 -0.5150854  2.71122936 0.4048951
# ALRU-ABPR -1.5137068 -3.3384178  0.31100413 0.1774172
# PSME-ABPR  1.2581284 -0.4452694  2.96152616 0.3035267
# TABR-ABPR  1.7344853 -0.0425305  3.51150118 0.0609028
# TSHE-ABPR  1.5843038 -0.1190939  3.28770161 0.0874328
# PSME-ALRU  2.7718352  1.0574913  4.48617906 0.0000488
# TABR-ALRU  3.2481922  1.4606810  5.03570333 0.0000027
# TSHE-ALRU  3.0980107  1.3836668  4.81235452 0.0000031
# TABR-PSME  0.4763570 -1.1871304  2.13984430 0.9794653
# TSHE-PSME  0.3261755 -1.2584282  1.91077914 0.9964747
# TSHE-TABR -0.1501815 -1.8136689  1.51330584 0.9999695

# significant differences are 

# ALRU-ABAM = 0.043 *
# ALRU-ABGR = 0.012 *
# PSME-ALRU = 0.00005 ***
# TABR-ALRU = 0.000003 ***
# TSHE-ALRU = 0.000003 ***

# ALRU is quite different from many of them 


# Nicer boxplot 
distances_EM2 <- data.frame(
  Host = betadisper.host$group,
  DistanceToCentroid = betadisper.host$distances
)

centroid_plot_EM2 <- ggplot(distances_EM2, aes(x = Host, y = DistanceToCentroid, fill = Host)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Distance to Centroid") +
  scale_fill_manual(values=all_hosts, 
                      name="Host",
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

centroid_plot_EM2


permanova.host <- vegan::adonis2(aitchison_EM ~ Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.host

# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# # 

# vegan::adonis2(formula = aitchison_EM ~ Host_ID, data = sample_data, permutations = 999, method = "euclidean")
#            Df SumOfSqs      R2      F Pr(>F)    
# Model      6    691.5 0.02779 1.7343  0.001 ***
# Residual 364  24188.5 0.97221                  
# Total    370  24879.9 1.00000 


# adonis result is significant but could be impacted by the heterogeneity in group variances observed 
# in the betadispersion 


########

# INTERPRETATION 

# There are significant differences in the beta dispersion of communities between host tree taxa, but not sites 
# The distance to centroid values for the tree communities are related to some environmental variables, but 
# simple linear relationships are an over-simplification here. 

## -- END -- ## 
