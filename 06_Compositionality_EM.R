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

clr_EM <- decostand(asv, method = "clr", pseudocount = 1e-06)

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
# Scaling because these CLR values are not scaled 
pca_clr_EM = prcomp(clr_sites_EM[7:1633], center = T, scale = T)

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
# Groups      3   3600 1199.95 5.2293    999  0.002 **
#   Residuals 367  84214  229.46                     
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
#                      diff         lwr       upr     p adj
# Northern-Andrews   4.630198  -0.9406522 10.201048 0.1409528
# Southern-Andrews  -2.351687  -8.0623340  3.358960 0.7123834
# WFDP-Andrews      -3.452371  -9.2333425  2.328601 0.4139034
# Southern-Northern -6.981885 -12.7063083 -1.257461 0.0095943
# WFDP-Northern     -8.082569 -13.8771496 -2.287988 0.0020446
# WFDP-Southern     -1.100684  -7.0297890  4.828421 0.9636782

# sig diffs between sites 

# southern - northern = 0.0096 *
# WFDP - northern = 0.002 * 


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
# Df SumOfSqs      R2      F Pr(>F)
# Model      3    32952 0.02542 3.1913  0.001 ***
#   Residual 367  1263152 0.97458
# Total    370  1296104 1.00000
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# adonis result is significant, and there was no heterogeneity of variances for the sites, so 
# this can be reported without concern 


permanova.both <- vegan::adonis2(aitchison_EM ~ Site * Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.both

# explained 11.8% 


permanova.both2 <- vegan::adonis2(aitchison_EM ~ Site + Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.both2

# explained 5.5%


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
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      6   9753 1625.44 7.7822    999  0.001 ***
#   Residuals 364  76027  208.87    
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
# ABGR-ABAM   2.1822157  -5.8168330 10.1812644 0.9839992
# ABPR-ABAM  -2.9693112 -11.4193589  5.4807365 0.9438988
# ALRU-ABAM  -9.1229818 -17.6277229 -0.6182407 0.0264743
# PSME-ABAM   5.2258020  -2.6305013 13.0821053 0.4342503
# TABR-ABAM   8.5380047   0.2873962 16.7886131 0.0371193
# TSHE-ABAM   4.3586037  -3.4976996 12.2149070 0.6531682
# ABPR-ABGR  -5.1515269 -13.7646131  3.4615592 0.5669858
# ALRU-ABGR -11.3051975 -19.9719483 -2.6384468 0.0024753
# PSME-ABGR   3.0435862  -4.9878177 11.0749901 0.9206193
# TABR-ABGR   6.3557890  -2.0617217 14.7732997 0.2774927
# TSHE-ABGR   2.1763880  -5.8550159 10.2077919 0.9845469
# ALRU-ABPR  -6.1536706 -15.2383332  2.9309919 0.4111236
# PSME-ABPR   8.1951132  -0.2855692 16.6757956 0.0658770
# TABR-ABPR  11.5073159   2.6601123 20.3545195 0.0025824
# TSHE-ABPR   7.3279149  -1.1527675 15.8085973 0.1409939
# PSME-ALRU  14.3487838   5.8136043 22.8839633 0.0000198
# TABR-ALRU  17.6609865   8.7615299 26.5604431 0.0000002
# TSHE-ALRU  13.4815855   4.9464060 22.0167650 0.0000813
# TABR-PSME   3.3122027  -4.9697782 11.5941836 0.8992725
# TSHE-PSME  -0.8671983  -8.7564422  7.0220457 0.9999030
# TSHE-TABR  -4.1794010 -12.4613819  4.1025799 0.7470465

# significant differences are 

# ALRU-ABAM = 0.026 *
# TABR-ABAM = 0.037 *
# ALRU-ABGR = 0.002 *
# TABR-ABPR = 0.002 * 
# PSME-ALRU = 0.00002 ***
# TABR-ALRU = 0.000002 ***
# TSHE-ALRU = 0.000008 ***

# ALRU is quite different from many of them, and TABR is different from a couple 


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
# Model      6    38664 0.02983 1.8654  0.001 ***
# Residual 364  1257440 0.97017                  
# Total    370  1296104 1.00000    


# adonis result is significant but could be impacted by the heterogeneity in group variances observed 
# in the betadispersion 


########

# INTERPRETATION 

# There are significant differences in the beta dispersion of communities between host tree taxa and sites 

## -- END -- ## 
