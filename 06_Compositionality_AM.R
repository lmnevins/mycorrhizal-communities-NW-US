# -----------------------------------------------------------------------------#
# Transforming AM phyloseq data and doing compositional analyses
# Original Author: L. McKinley Nevins 
# April 21, 2025
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
#                     cowplot v 1.1.3
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
library(cowplot); packageVersion("cowplot")


#################################################################################
#                               Main workflow                                   #
#  Transform data to account for compositionality, then perform appropriate     #
#  analyses that allow for calculating distance values. 
#                                                                               #
#################################################################################

## If just running analyses, skip to STEP 2 

#################################################################################

################ --
# (1) DATA PREP
################ --

## GOAL: Create cleaned and transformed phyloseq object for all downstream analyses 


# Load phyloseq object produced from the subset of the funguild community 
ps_AM <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_final_2025.RDS")

# Trim down just to include taxa that had functional assignments made 

# Load in taxon table for AM community
tax <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_subset_tax_2025.csv")

# pull out the genus and family levels, as well as the original OTU's. This will be necessary 
# to match them up with the site x environment matrix 

tax <- dplyr::select(tax, X, Family, Genus)

tax$Family <- as.factor(tax$Family)

summary(tax)


# Load in table of Families with functional classifications 

exp <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/AM_functions_all.csv")

exp$Family <- as.factor(exp$Family)

summary(exp)


# Merge datasets together 
matrix <- merge(tax, exp, by = "Family", all.y = TRUE)

#get column name set for name merging down below
matrix$ASV <- matrix$X

# the number is correct - 363 ASVs with genus and functional assignments 

#start names file with the original ASVs 
ASV <- phyloseq::taxa_names(ps_AM)
ASV_long <- as.data.frame(ASV)

#change OTU names to something nicer to work with
taxa_names(ps_AM)
n_seqs <- seq(ntaxa(ps_AM))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_AM) <- paste("ASV", formatC(n_seqs, 
                                              width = len_n_seqs, 
                                              flag = "0"), sep = "_")
taxa_names(ps_AM)

# get shortened names 
ASV2 <- taxa_names(ps_AM) 
ASV_short <- as.data.frame(ASV2)

# join two dataframes
all_names <- cbind(ASV_long, ASV_short)


# save all_names file as key to the long and short OTU assignments 
write.csv(all_names, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/all_ASV_names_AM.csv")


#merge taxa_names with the traits file to get the updated ASV names 
matrix <- merge(matrix, all_names, by = "ASV")

# Prune the final phyloseq object to just contain the taxa that have functional assignments 
taxa_to_keep <- matrix$ASV2

# Prune taxa from the phyloseq object that are NOT in the taxa_to_keep list 
ps_AM_trimmed_funcs <- prune_taxa(taxa_to_keep, ps_AM)

# pull out values 
AM_final <- otu_table(ps_AM_trimmed_funcs) %>% as("matrix") %>% as.data.frame()

# load this subset of taxa back into a phyloseq object 
ps_AM_final <- phyloseq(tax_table(tax_table(ps_AM)),
                        otu_table(otu_table(AM_final, taxa_are_rows=FALSE)),
                        sample_data(ps_AM))


# Getting to know your phyloseq data #### -- 

# number of taxa - 363
ntaxa(ps_AM_final)

# number of samples - 131
nsamples(ps_AM_final)

# sample names
sample_names(ps_AM_final)
rank_names(ps_AM_final)
# taxa names
taxa_names(ps_AM_final)


# ASV table
otu_table(ps_AM_final) %>% View()

# how many sequences observed in each sample?
seq_counts <- otu_table(ps_AM_final) %>% rowSums() %>% as.data.frame()


####### -- 
# Note: There are 16 trees with <10 ASV's
###### -- 


# Remove one tree, S-TABR-01, that doesn't have any reads 
ps_AM_final <- subset_samples(ps_AM_final, sample_names(ps_AM_final) != "S-TABR-01_FWD_filt.fastq.gz")

# Check again
otu_table(ps_AM_final) %>% View() # Gone, now 130 trees 


# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_AM_final) %>% colSums() %>% as.data.frame()


# how many different samples was each taxon found in?
asv <- otu_table(ps_AM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

#inspect tax_table 

AM_tax_table <- as.data.frame(tax_table(ps_AM_final))

#count number of classifications in each column to determine coverage to taxonomic levels 

colSums(!is.na(AM_tax_table))

#  Kingdom  Phylum    Class    Order   Family   Genus    Species 
#  363      363       363      363     363      296       197
#  100%     100%      100%     100%    100%     81.5%    54.3% 


# Load in sample data for the communities 
sample_metadata <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_enviro_all_2025.csv")

sample_metadata <- column_to_rownames(sample_metadata, "X")


############# -- 
# Save for untransformed final phyloseq 

saveRDS(ps_AM_final, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_func_subset_final_2025.RDS")

# This should be used for the differential abundance analyses later on 


############# -- 
# Save 'matrix' which now has functional info matched up to asvs and taxa
write.csv(matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_ASV_funcs_AM.csv")

############# -- 


####### more data exploration ######### -- 

#check rarefaction curve 
#using raw counts asv table 

otu.rare = otu_table(ps_AM_final)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

# we will use vegan rarecurve 
otu.rarecurve = rarecurve(otu.rare, step = 10000, label = F)

# very different sample depths 

## The analyses in the next section are using the raw data 
# will compare this with the transformed data that will be generated below to see if the 
# transformation has a significant effect on the composition 

#### 

# Pursuing a transformation because sequence data are inherently compositional, but standard rarefaction 
# is not repeatable and leads to the exclusion of a huge fraction of your data 

# Centered log-ratio transformation per https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full

# Also: Van den Boogaart, K. G., and Tolosana-Delgado, R. (2013). Analyzing Compositional Data with R, 
# London, UK: Springer.


## Using clr transformation with a small pseudocount to account for many zeros in the data 

clr_AM <- decostand(asv, method = "clr",  pseudocount = 1e-06)

clr_AM <- as.data.frame(clr_AM)

# load the clr data back into a phyloseq object 
ps_AM_clr <- phyloseq(tax_table(tax_table(ps_AM_final)),
                     otu_table(otu_table(clr_AM, taxa_are_rows=FALSE)),
                     sample_data(sample_metadata))


##################### -- 
# save phyloseq object with transformed data 
saveRDS(ps_AM_clr, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")

## This is the phyloseq object to be used in all downstream analyses 
##################### -- 

#################################################################################

## SKIP TO HERE ## 

######################################## --
# (2) ANALYSES OF COMMUNITY COMPOSITION
######################################## --

## GOAL: Analyze the taxonomic composition of the AM communities.  


# Load in ps_AM_clr dataframe 
ps_AM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")


# Pull out otu table to generate 'clr_AM' dataframe 
clr_AM <- otu_table(ps_AM_clr) %>% as.matrix() %>% as.data.frame()


## Beta-diversity can be calculated using Aitchison Distance, which is essentially the euclidian
# distance calculated between pairs of samples that have been transformed by CLR

# calculate Aitchison distance using dist() from base R 

aitchison_AM <- dist(clr_AM, method = "euclidean")

# save Aitchison Distance matrix 
save(aitchison_AM, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_AM_2025.Rdata")


## visualize ##

# PCA is appropriate for euclidean distances 

# get a few site-level columns to load back into the clr dataframe 

sites <- dplyr::select(sample_metadata, Sample_ID, Site, Host_ID, Field_ID, Location)

clr_sites_AM <- merge(sites, clr_AM, by = 'row.names')

#PCA of differences in composition for Sites 
pca_clr_AM = prcomp(clr_sites_AM[7:369], center = T, scale = T)

sd.pca_clr_AM = pca_clr_AM$sdev
loadings.pca_clr_AM = pca_clr_AM$rotation
names.pca_clr_AM = colnames(clr_sites_AM[7:369])
scores.pca_clr_AM = as.data.frame(pca_clr_AM$x)
scores.pca_clr_AM$Site = clr_sites_AM$Site
scores.pca_clr_AM$Host_ID = clr_sites_AM$Host_ID
summary(pca_clr_AM)



# PCA scores are 'scores.pca_clr_AM' with column for Sites

loadings.pca_clr_AM <- as.data.frame(loadings.pca_clr_AM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_clr_AM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# For just AM hosts 
AM_hosts <- c("#B93289FF", "#F48849FF", "#ffe24cFF")

sites <- c(15,16,17,18)

# Plot the Results by site and host 
PCA_plot_both <- ggplot(scores.pca_clr_AM, aes(x = PC1, y = PC2, color = Host_ID, shape = Site)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_colour_manual(values=AM_hosts, 
                      name="Host Tree Species",
                      breaks=c("ALRU", "TABR", "THPL"),
                      labels=c("ALRU", "TABR", "THPL")) +
  scale_shape_manual(values=sites, 
                     name="Site",
                     breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                     labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  labs(title = "PCA Biplot: Community Variation Across Sites",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
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

PCA_plot_both


# Plot the Results by site alone
PCA_plot_site <- ggplot(scores.pca_clr_AM, aes(x = PC1, y = PC2, color = Site)) +
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


## Try out a way to visualize differences in the communities along these two axes: 


### Calculate average PC1 score for each site, compare statistically  ####

PC1_site <- dplyr::select(scores.pca_clr_AM, PC1, Site, Host_ID)

PC1_site_summary <- PC1_site %>%
  group_by(Site) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1   = sd(PC1, na.rm = TRUE),
    n        = n(),
    se_PC1   = sd_PC1 / sqrt(n)
  )


# Test for significant differences between species 
aov_PC1_site <- aov(PC1 ~ Site, data = PC1_site)
summary(aov_PC1_site) # Significant

tuk_PC1_site <- TukeyHSD(aov_PC1_site)
tuk_PC1_site

# Significant differences between sites 
# Northern-Andrews p = 0.0000008
# Southern-Northern p = 0.0000010
# WFDP-Northern p = 0.0055585


# Visualize 

# Define preferred order
site_order <- c("Northern", "WFDP", "Andrews", "Southern")

# Apply to site data
PC1_site_summary$Site <- factor(PC1_site_summary$Site, levels = site_order)


PC1_site_plot <- ggplot(PC1_site_summary, aes(x = Site, y = mean_PC1, fill = Site)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values=palette, 
                    name="Site",
                    breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                    labels=c("Northern", "WFDP", "Andrews", "Southern")) + 
  labs(title = "", x = "", y = "PCA Axis 1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC1_site_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 


### Calculate average PC2 score for each site, compare statistically 

PC2_site <- dplyr::select(scores.pca_clr_AM, PC2, Site, Host_ID)

PC2_site_summary <- PC2_site %>%
  group_by(Site) %>%
  summarise(
    mean_PC2 = mean(PC2, na.rm = TRUE),
    sd_PC2   = sd(PC2, na.rm = TRUE),
    n        = n(),
    se_PC2   = sd_PC2 / sqrt(n)
  )

# Test for significant differences between sites
aov_PC2_site <- aov(PC2 ~ Site, data = PC2_site)
summary(aov_PC2_site) # Significant

tuk_PC2_site <- TukeyHSD(aov_PC2_site)
tuk_PC2_site

# Significant differences between sites
# Southern-Andrews p = 0.0227669
# Southern-Northern p = 0.0455868


# Visualize 

# Apply order to site data
PC2_site_summary$Site <- factor(PC2_site_summary$Site, levels = site_order)


PC2_site_plot <- ggplot(PC2_site_summary, aes(x = Site, y = mean_PC2, fill = Site)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values=palette, 
                    name="Site",
                    breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                    labels=c("Northern", "WFDP", "Andrews", "Southern")) + 
  labs(title = "", x = "", y = "PCA Axis 2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC2_site_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 

######################### 

#### -- 

# Plot the Results by host alone
PCA_plot_host <- ggplot(scores.pca_clr_AM, aes(x = PC1, y = PC2, color = Host_ID)) +
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

PCA_plot_host



### Calculate average PC1 score for each host, compare statistically   ######

PC1_host <- dplyr::select(scores.pca_clr_AM, PC1, Site, Host_ID)

PC1_host_summary <- PC1_host %>%
  group_by(Host_ID) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1   = sd(PC1, na.rm = TRUE),
    n        = n(),
    se_PC1   = sd_PC1 / sqrt(n)
  )


# Test for significant differences between hosts 
aov_PC1_host <- aov(PC1 ~ Host_ID, data = PC1_host)
summary(aov_PC1_host) # NOT Significant

tuk_PC1_host <- TukeyHSD(aov_PC1_host)
tuk_PC1_host

# NO Significant differences between hosts



# Visualize 

PC1_host_plot <- ggplot(PC1_host_summary, aes(x = Host_ID, y = mean_PC1, fill = Host_ID)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values=AM_hosts, 
                    name="Host Tree Species",
                    breaks=c("ALRU", "TABR", "THPL"),
                    labels=c("ALRU", "TABR", "THPL")) + 
  labs(title = "", x = "", y = "PCA Axis 1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC1_host_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 


### Calculate average PC2 score for each host, compare statistically 

PC2_host <- dplyr::select(scores.pca_clr_AM, PC2, Site, Host_ID)

PC2_host_summary <- PC2_host %>%
  group_by(Host_ID) %>%
  summarise(
    mean_PC2 = mean(PC2, na.rm = TRUE),
    sd_PC2   = sd(PC2, na.rm = TRUE),
    n        = n(),
    se_PC2   = sd_PC2 / sqrt(n)
  )

# Test for significant differences between hosts
aov_PC2_host <- aov(PC2 ~ Host_ID, data = PC2_host)
summary(aov_PC2_host) # Significant

tuk_PC2_host <- TukeyHSD(aov_PC2_host)
tuk_PC2_host

# Significant difference between hosts
# TABR-ALRU p = 0.0227506


# Visualize 

PC2_host_plot <- ggplot(PC2_host_summary, aes(x = Host_ID, y = mean_PC2, fill = Host_ID)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values=AM_hosts, 
                    name="Host Tree Species",
                    breaks=c("ALRU", "TABR", "THPL"),
                    labels=c("ALRU", "TABR", "THPL")) + 
  labs(title = "", x = "", y = "PCA Axis 2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC2_host_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 

#### -- 



#####assess multivariate homogeneity of sites###########
#betadisper function in vegan 

#first object needs to be a dist object of Aitchison distances 
#second object is the groups of interest, as a vector
betadisper.site <- vegan::betadisper(aitchison_AM, sample_metadata$Site, type = "median", sqrt.dist = FALSE)

betadisper.site

# Permutation test for homogeneity of multivariate dispersions
# p-value > 0.05, there is no evidence the groups have different dispersions
# p-value < 0.05, there is evidence the groups have different dispersions 

# This is permutational, so it will give a different p-value each time 
permutest(betadisper.site)
#groups dispersions are different (not homogenous) p = 0.003

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
#                        diff        lwr       upr     p adj
# Northern-Andrews   -8.495812 -24.354237  7.362613 0.5049208
# Southern-Andrews    2.174041 -14.979765 19.327847 0.9875469
# WFDP-Andrews      -21.897533 -38.668678 -5.126388 0.0049518
# Southern-Northern  10.669853  -7.238962 28.578668 0.4103049
# WFDP-Northern     -13.401721 -30.944352  4.140910 0.1975852
# WFDP-Southern     -24.071574 -42.793415 -5.349733 0.0058530

# differences between 
# WFDP-Andrews = 0.005 *
# WFDP-Southern = 0.006 * 


# Nicer boxplot 
distances_AM <- data.frame(
  Site = betadisper.site$group,
  DistanceToCentroid = betadisper.site$distances
)

# Define preferred order
site_order <- c("Northern", "WFDP", "Andrews", "Southern")

# Apply to distances data
distances_AM$Site <- factor(distances_AM$Site, levels = site_order)

centroid_plot_AM <- ggplot(distances_AM, aes(x = Site, y = DistanceToCentroid, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal() +
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



permanova.site <- vegan::adonis2(aitchison_AM ~ Site, data = sample_metadata, method = "euclidean", permutations = 999)
permanova.site

# vegan::adonis2(formula = aitchison_AM ~ Site, data = sample_metadata, permutations = 999, method = "euclidean")
# Df SumOfSqs      R2      F Pr(>F)    
#   Model      3    62951 0.09633 4.4773  0.001 ***
#Residual 126   590521 0.90367                  
# Total    129   653472 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# adonis result is significant but could be impacted by the heterogeneity in group variances observed 
# in the betadispersion 


permanova.both <- vegan::adonis2(aitchison_AM ~ Site * Host_ID, data = sample_metadata, method = "euclidean", permutations = 999)
permanova.both

# explained 21.6%


permanova.both2 <- vegan::adonis2(aitchison_AM ~ Site + Host_ID, data = sample_metadata, method = "euclidean", permutations = 999)
permanova.both2

# explained 14.3%


##############

#check betadispersion between host tree taxa 
betadisper.host <- betadisper(aitchison_AM, sample_metadata$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper.host

permutest(betadisper.host)
#groups are different p = 0.001

#if the dispersion is different between groups, then examine
plot(betadisper.host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#B93289FF", "#F48849FF", "#FEBC2AFF"))

boxplot(betadisper.host)
mod.HSD.host <- TukeyHSD(betadisper.host)
mod.HSD.host
plot(mod.HSD.host)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#              diff        lwr       upr      p adj
# TABR-ALRU -26.74842 -39.720975 -13.775865 0.0000089
# THPL-ALRU -10.85403 -23.481131   1.773063 0.1072439
# THPL-TABR  15.89439   3.347832  28.440940 0.0089417

# TABR - ALRU 0.00000008 ***
# TABR - THPL 0.0089 *

# Nicer boxplot 
distances_AM2 <- data.frame(
  Host = betadisper.host$group,
  DistanceToCentroid = betadisper.host$distances
)

centroid_plot_AM2 <- ggplot(distances_AM2, aes(x = Host, y = DistanceToCentroid, fill = Host)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal() +
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


permanova.host <- vegan::adonis2(aitchison_AM ~ Host_ID, data = sample_metadata, method = "euclidean", permutations = 999)
permanova.host


########

# INTERPRETATION 

# There are significant differences in the beta dispersion of communities between sites and host tree taxa 
# The distance to centroid values for the tree communities are related to some environmental variables, but 
# simple linear relationships are an over-simplification here. 


## -- END -- ## 

