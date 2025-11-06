# -----------------------------------------------------------------------------#
# Transforming EM phyloseq data and doing compositional analyses
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
#  Transform data to account for compositionality, then perform appropriate     #
#  analyses that allow for calculating distance values. 
#                                                                               #
#################################################################################

# Load phyloseq object produced from the subset of the funguild community 
ps_EM <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_final.RDS")

# Trim down just to include taxa that had functional assignments made 

# Load in taxon table for EM community
tax <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_subset_tax.csv")

# pull out the genus and family levels, as well as the original ASV's. 
tax <- dplyr::select(tax, X, Family, Genus)

tax$Genus <- as.factor(tax$Genus)

summary(tax)


# Load in table of Genera with functional classifications 
exp <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/em_functions_all.csv")

exp$Genus <- gsub("g__", "", exp$Genus)

exp$Genus <- as.factor(exp$Genus)

summary(exp)

# Merge datasets together 
matrix <- merge(tax, exp, by = "Genus", all.y = TRUE)

#get column name set for name merging down below
matrix$OTU <- matrix$X

# the number is correct - 1,941 OTUs with genus and functional assignments 

#start names file with the original OTUs 
OTU <- phyloseq::taxa_names(ps_EM)
OTU_long <- as.data.frame(OTU)

#change OTU names to something nicer to work with
taxa_names(ps_EM)
n_seqs <- seq(ntaxa(ps_EM))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_EM) <- paste("OTU", formatC(n_seqs, 
                                              width = len_n_seqs, 
                                              flag = "0"), sep = "_")
taxa_names(ps_EM)

# get shortened names 
OTU2 <- taxa_names(ps_EM) 
OTU_short <- as.data.frame(OTU2)

# join two dataframes
all_names <- cbind(OTU_long, OTU_short)

# save all_names file as key to the long and short OTU assignments 
write.csv(all_names, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/all_OTU_names_EM.csv")


#merge taxa_names with the traits file to get the updated OTU names 
matrix <- merge(matrix, all_names, by = "OTU")


# Prune the final phyloseq object to just contain the taxa that have functional assignments 
taxa_to_keep <- matrix$OTU2

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

# number of taxa - 1,941
ntaxa(ps_EM_final)

# number of samples - 425
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
# Note: There are 19 trees with <10 ASV's -> for due diligence should run the analyses without 
# them and see if it makes a difference. For now will retain. 
######

# Remove one tree, S-ALRU-05, that doesn't have any reads 
ps_EM_final <- subset_samples(ps_EM_final, sample_names(ps_EM_final) != "S-ALRU-05_FWD")

# Check again
otu_table(ps_EM_final) %>% View() # Gone, now 424 trees 


# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_EM_final) %>% colSums()


# how many different samples was each taxon found in?
asv <- otu_table(ps_EM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

#inspect tax_table 

EM_tax_table <- as.data.frame(tax_table(ps_EM_final))

#count number of classifications in each column to determine coverage to taxonomic levels 

colSums(!is.na(EM_tax_table))

#  Kingdom  Phylum    Class    Order   Family   Genus    Species 
#  1941     1941      1941     1941    1941     1941     905
#  100%     100%      100%     100%    100%     100%     46.6% 


# Load in sample data for the communities 
sample_metadata <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_all.csv")

sample_metadata <- column_to_rownames(sample_metadata, "X")


#############
# Save for untransformed final phyloseq 

saveRDS(ps_EM_final, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_func_subset_final.RDS")

# This should be used for the differential abundance analyses later on 

#############

# Save 'matrix' which now has functional info matched up to OTU's and taxa
write.csv(matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_OTU_funcs_EM.csv")


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
saveRDS(ps_EM_clr, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final.RDS")

## This is the phyloseq object to be used in all downstream analyses 
#####################

#####################Beta-diversity##########################################

## Beta-diversity can be calculated using Aitchison Distance, which is essentially the euclidean
# distance calculated between pairs of samples that have been transformed by CLR

# calculate Aitchison distance using dist() from base R 

aitchison_EM <- dist(clr_EM, method = "euclidean")

# save Aitchison Distance matrix 
save(aitchison_EM, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_2025.Rdata")


## visualize ##

# PCA is appropriate for euclidean distances 

# get a few site-level columns to load back into the clr dataframe 

sample_data <-  data.frame(sample_data(ps_EM_clr))

sites <- dplyr::select(sample_data, Sample_ID, Site, Host_ID, Field_ID, Location)

clr_sites_EM <- merge(sites, clr_EM, by = 'row.names')

#PCA of differences in composition for Sites 
pca_clr_EM = prcomp(clr_sites_EM[7:1947], center = T, scale = F)

sd.pca_clr_EM = pca_clr_EM$sdev
loadings.pca_clr_EM = pca_clr_EM$rotation
names.pca_clr_EM = colnames(clr_sites_EM[7:1947])
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
# ABAM        ABGR        ABPR          ALRU         PSME         TABR         THPL       TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66")

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
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE")) +
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
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE")) +
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
# Groups      3   63.2 21.0721 2.2334    999  0.081 .
# Residuals 420 3962.7  9.4351                       
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
#                             diff                     lwr                     upr                   p adj
# Northern-Andrews   0.59326323424309546084 -0.46321462583038353067 1.649741094316574452350 0.469877967663406836962
# Southern-Andrews  -0.17059213228115321925 -1.25902797535443067289 0.917843710792124234388 0.977640988003998123723
# WFDP-Andrews      -0.46137057741935638688 -1.54688979709627005832 0.624148642257557284552 0.692011701597894823834
# Southern-Northern -0.76385536652424868009 -1.85910626201394224566 0.331395528965444885472 0.275282230178262810050
# WFDP-Northern     -1.05463381166245184772 -2.14698628036122984852 0.037718657036326153076 0.062876738357810046942
# WFDP-Southern     -0.29077844513820316763 -1.41406871551960056088 0.832511825243194225621 0.909207565909160075890

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
  labs(title = "Beta Dispersion by Site",
       x = "Site",
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
#           Df   SumOfSqs   R2     F  Pr(>F)    
# Model      3      469 0.0133 1.892  0.001 ***
# Residual 421    34780 0.9867                 
# Total    424    35249 1.0000                  
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

# When including the interaction of tree host and site, 9.6% of the variaiton was explained. Envrionmental and host identity explain only 
# 9.6% of the variation. Only enviro explained 1%, ... just compare 


# One thing for the discussion is that there's low variation explained. Maybe fungi are not that sensitive to changes in 
# the environment. 


# Not finding a lot of variation explained is interesting. 


##################

# Take distance to centroid values and regress with key environmental variables to assess if there are 
# relationships with environmental variation across the sites

# need to make dataset that has distance to centroid and the environmental variables 

# Have distance to centroid and environmental data for each individual tree 
distances_EM <- tibble::rownames_to_column(distances_EM, "Tree")


# Sample metadata has all of the environmental variables of interest 
EM_env <- rownames_to_column(sample_metadata, "Tree")

# pick some specific environmental variables of interest to compare 
EM_env2 <- dplyr::select(EM_env, Tree, elev, mean_precip_mm, mean_summer_precip_mm, MAT, pct_N, pct_C, Sand, 
                         ph, org_matter, EC, avg_July_SPEI, count_mod_dry, count_sev_dry, apr1_SWE)

# combine with the distance to centroid summary table
centroid_enviro_EM <- merge(distances_EM, EM_env2, by = "Tree")


### Testing these using the individual tree data

## Elevation
elev2 <- ggplot(centroid_enviro_EM, aes(x = elev, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
elev2

# test relationships 
lm_elev2 <- lm(DistanceToCentroid ~ elev, data = centroid_enviro_EM)
summary(lm_elev2) # Not significant 


## Mean annual precip
MAP2 <- ggplot(centroid_enviro_EM, aes(x = mean_precip_mm, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MAP2

# test relationships 
lm_MAP2 <- lm(DistanceToCentroid ~ mean_precip_mm, data = centroid_enviro_EM)
summary(lm_MAP2) # Not significant 

## Mean summer precip
MSP2 <- ggplot(centroid_enviro_EM, aes(x = mean_summer_precip_mm, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MSP2

# test relationships 
lm_MSP2 <- lm(DistanceToCentroid ~ mean_summer_precip_mm, data = centroid_enviro_EM)
summary(lm_MSP2) # SIGNIFICANT 


# Adjusted R-squared:  0.013713717639779022 , p-value: 0.0090244468063364527


## MAT
MAT2 <- ggplot(centroid_enviro_EM, aes(x = MAT, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MAT2

# test relationships 
lm_MAT2 <- lm(DistanceToCentroid ~ MAT, data = centroid_enviro_EM)
summary(lm_MAT2) # not significant 


## Soil percent N
pctN2 <- ggplot(centroid_enviro_EM, aes(x = pct_N, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
pctN2

# test relationships 
lm_pctN2 <- lm(DistanceToCentroid ~ pct_N, data = centroid_enviro_EM)
summary(lm_pctN2) # not significant 


## Soil percent C
pctC2 <- ggplot(centroid_enviro_EM, aes(x = pct_C, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
pctC2

# test relationships 
lm_pctC2 <- lm(DistanceToCentroid ~ pct_C, data = centroid_enviro_EM)
summary(lm_pctC2) # not significant 


## Soil Organic Matter
SOM2 <- ggplot(centroid_enviro_EM, aes(x = org_matter, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SOM2

# test relationships 
lm_SOM2 <- lm(DistanceToCentroid ~ org_matter, data = centroid_enviro_EM)
summary(lm_SOM2) # not significant 


## Soil sand fraction
sand2 <- ggplot(centroid_enviro_EM, aes(x = Sand, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
sand2

# test relationships 
lm_sand2 <- lm(DistanceToCentroid ~ Sand, data = centroid_enviro_EM)
summary(lm_sand2) # SIGNIFICANT 

# Adjusted R-squared:  0.04909630354015726 ,  p-value: 0.0000024349251781372014


## pH
ph2 <- ggplot(centroid_enviro_EM, aes(x = ph, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ph2

# test relationships 
lm_ph2 <- lm(DistanceToCentroid ~ ph, data = centroid_enviro_EM)
summary(lm_ph2) # not signfiicant 


## Cation exchange capacity (EC)
EC2 <- ggplot(centroid_enviro_EM, aes(x = EC, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
EC2

# test relationships 
lm_EC2 <- lm(DistanceToCentroid ~ EC, data = centroid_enviro_EM)
summary(lm_EC2) # not significant 


## Average July SPEI (avg_July_SPEI)
SPEI2 <- ggplot(centroid_enviro_EM, aes(x = avg_July_SPEI, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SPEI2

# test relationships 
lm_SPEI2 <- lm(DistanceToCentroid ~ avg_July_SPEI, data = centroid_enviro_EM)
summary(lm_SPEI2) # not significant 


## Count of Moderate Dry Months (count_mod_dry)
MOD2 <- ggplot(centroid_enviro_EM, aes(x = count_mod_dry, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MOD2

# test relationships 
lm_MOD2 <- lm(DistanceToCentroid ~ count_mod_dry, data = centroid_enviro_EM)
summary(lm_MOD2) # SIGNIFICANT 

# Adjusted R-squared:  0.024987687527389557 ,  p-value: 0.00063727795318942832


## Count of Severe Dry Months (count_sev_dry)
SEV2 <- ggplot(centroid_enviro_EM, aes(x = count_sev_dry, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SEV2

# test relationships 
lm_SEV2 <- lm(DistanceToCentroid ~ count_sev_dry, data = centroid_enviro_EM)
summary(lm_SEV2) # SIGNIFICANT 

# Adjusted R-squared:  0.015578269645406762 ,  p-value: 0.0057866272483853729


# Are moderate and severe counts related? 
mod_sev <- ggplot(centroid_enviro_EM, aes(x = count_sev_dry, y = count_mod_dry)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
mod_sev


# test relationships 
lm_mod_sev <- lm(count_mod_dry ~ count_sev_dry, data = centroid_enviro_EM)
summary(lm_mod_sev) 
# They are related but it's not super strong 


## Average April 1 Snow Water Equivalent (apr1_SWE)
SWE2 <- ggplot(centroid_enviro_EM, aes(x = apr1_SWE, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SWE2

# test relationships 
lm_SWE2 <- lm(DistanceToCentroid ~ apr1_SWE, data = centroid_enviro_EM)
summary(lm_SWE2) # not significant 



### PLOT SIGNIFICANT RELATIONSHIPS ###

# Factors related to soil and drought 

#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# removed the legends from these so they will plot better as panels. Can add back in using theme(legend.position = "right")


# Mean Annual Precipitation
MSP_plot <- ggplot(centroid_enviro_EM, aes(x = mean_summer_precip_mm, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40") +
  theme_minimal() +
  labs(x = "Mean Summer Precipitation (mm)", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

MSP_plot


# Soil Sand Fraction 
soil_Sand_EM_plot <- ggplot(centroid_enviro_EM, aes(x = Sand, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Soil Sand (%)", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

soil_Sand_EM_plot


# Count of Moderate Dry Months 
MOD_dry_plot <- ggplot(centroid_enviro_EM, aes(x = count_mod_dry, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Count of Moderately Dry Months", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

MOD_dry_plot


# Count of Severe Dry Months 
SEV_dry_plot <- ggplot(centroid_enviro_EM, aes(x = count_sev_dry, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Count of Severely Dry Months", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

SEV_dry_plot


# Could do a nice multi-panel figure with these 
cowplot::plot_grid(MSP_plot, soil_Sand_EM_plot, MOD_dry_plot, SEV_dry_plot, nrow = 2, ncol = 2)


##############

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
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      7  330.8  47.263 5.5804    999  0.001 ***
#   Residuals 416 3523.3   8.470                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#if the dispersion is different between groups, then examine
plot(betadisper.host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66"))

boxplot(betadisper.host)
mod.HSD.host <- TukeyHSD(betadisper.host)
mod.HSD.host
plot(mod.HSD.host)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
#                     diff                      lwr                    upr                      p adj
# ABGR-ABAM  0.24009000765485133400 -1.414969624121183544574 1.89514963943088621257 0.999850164022045628797741
# ABPR-ABAM -0.21088470085974719836 -1.959259199033201115370 1.53748979731370671864 0.999956972022758572293810
# ALRU-ABAM -1.42023727182281511716 -3.179928227929698980603 0.33945368428406874628 0.216372438285674428826155
# PSME-ABAM  1.04368502706741406172 -0.581839569762411024101 2.66920962389723914754 0.513125793913917327415675
# TABR-ABAM  1.56342326456266089707 -0.143685843785576583542 3.27053237291089837768 0.100275713314807446430166
# THPL-ABAM  0.75107664882099811621 -0.920196199069523323999 2.42234949671151955641 0.870816059834571998266028
# TSHE-ABAM  1.43079986915775503320 -0.194724727672070052620 3.05632446598758011902 0.131171454008395249601904
# ABPR-ABGR -0.45097470851459853236 -2.233083013665135840853 1.33113359663593877613 0.994475427855424976009147
# ALRU-ABGR -1.66032727947766645116 -3.453539170590259388405 0.13288461163492670813 0.092683406013497116049393
# PSME-ABGR  0.80359501941256272772 -0.858159126158318352751 2.46534916498344358615 0.821391632048436459712093
# TABR-ABGR  1.32333325690780956307 -0.418309109718806304556 3.06497562353442543071 0.287991557355979366228382
# THPL-ABGR  0.51098664116614678221 -1.195544799364422150489 2.21751808169671571491 0.984826351852879700032872
# TSHE-ABGR  1.19070986150290369920 -0.471044284067977381270 2.85246400707378455763 0.364254087836324069726857
# ALRU-ABPR -1.20935257096306791880 -3.089033365022595400973 0.67032822309645956338 0.510365364689852052393348
# PSME-ABPR  1.25456972792716126008 -0.500143308854255508322 3.00928276470857802849 0.367184750823381644835308
# TABR-ABPR  1.77430796542240809544 -0.056240896082635494224 3.60485682692745168509 0.065091663663475074486087
# THPL-ABPR  0.96196134968074531457 -0.835214367562282289370 2.75913706692377314056 0.731488717268277133243259
# TSHE-ABPR  1.64168457001750223156 -0.113028466763914536841 3.39639760679891899997 0.085817099904708427082767
# PSME-ALRU  2.46392229889022917888  0.697933421068020587441 4.22991117671243799236 0.000689481828908267146971
# TABR-ALRU  2.98366053638547601423  1.142300167120813725319 4.82502090565013830314 0.000031534084967477227224
# THPL-ALRU  2.17131392064381323337  0.363127139973997259403 3.97950070131362920733 0.006894588230373699389020
# TSHE-ALRU  2.85103714098057015036  1.085048263158361558922 4.61702601880277896385 0.000034393494188877937745
# TABR-PSME  0.51973823749524683535 -1.193862055278099854050 2.23333853026859330271 0.983647470010056590439262
# THPL-PSME -0.29260837824641594551 -1.970511053517981947536 1.38529429702515005651 0.999489719953239097449682
# TSHE-PSME  0.38711484209034097148 -1.245225401390211539265 2.01945508557089326018 0.996307142157155212203179
# THPL-TABR -0.81234661574166278086 -2.569403439122804222450 0.94471020763947888277 0.853155847853298676675138
# TSHE-TABR -0.13262339540490586387 -1.846223688178252553271 1.58097689736844082553 0.999997957060200781675974
# TSHE-THPL  0.67972322033675691699 -0.998179454934809085032 2.35762589560832314106 0.921360265450602122783152

# significant differences are 

# PSME-ALRU
# TABR-ALRU
# THPL-ALRU
# TSHE-ALRU

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
  labs(title = "Beta Dispersion by Host",
       x = "Host",
       y = "Distance to Centroid") +
  scale_fill_manual(values=all_hosts, 
                      name="Host",
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE")) +
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
#           Df              SumOfSqs                    R2                   F Pr(>F)    
# Model      7   818.660247541602871 0.0275682523534271327 1.68479000000000001  0.001 ***
# Residual 416 28877.101277212157584 0.9724317476465739185                               
# Total    423 29695.761524753728736 1.0000000000000000000                               
# ---


# adonis result is significant but could be impacted by the heterogeneity in group variances observed 
# in the betadispersion 


########

# INTERPRETATION 

# There are significant differences in the beta dispersion of communities between host tree taxa, but not sites 
# The distance to centroid values for the tree communities are related to some environmental variables, but 
# simple linear relationships are an over-simplification here. 


# Save centroid_enviro_EM dataframe to use in future analyses 

write.csv(centroid_enviro_EM, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_centroid_distance_enviro.csv")



