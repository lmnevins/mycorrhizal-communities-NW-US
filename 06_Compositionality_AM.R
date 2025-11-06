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
matrix$OTU <- matrix$X

# the number is correct - 363 OTUs with genus and functional assignments 

#start names file with the original OTUs 
OTU <- phyloseq::taxa_names(ps_AM)
OTU_long <- as.data.frame(OTU)

#change OTU names to something nicer to work with
taxa_names(ps_AM)
n_seqs <- seq(ntaxa(ps_AM))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_AM) <- paste("OTU", formatC(n_seqs, 
                                              width = len_n_seqs, 
                                              flag = "0"), sep = "_")
taxa_names(ps_AM)

# get shortened names 
OTU2 <- taxa_names(ps_AM) 
OTU_short <- as.data.frame(OTU2)

# join two dataframes
all_names <- cbind(OTU_long, OTU_short)


# save all_names file as key to the long and short OTU assignments 
write.csv(all_names, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/all_OTU_names_AM.csv")


#merge taxa_names with the traits file to get the updated OTU names 
matrix <- merge(matrix, all_names, by = "OTU")

# Prune the final phyloseq object to just contain the taxa that have functional assignments 
taxa_to_keep <- matrix$OTU2

# Prune taxa from the phyloseq object that are NOT in the taxa_to_keep list 
ps_AM_trimmed_funcs <- prune_taxa(taxa_to_keep, ps_AM)

# pull out values 
AM_final <- otu_table(ps_AM_trimmed_funcs) %>% as("matrix") %>% as.data.frame()

# load this subset of taxa back into a phyloseq object 
ps_AM_final <- phyloseq(tax_table(tax_table(ps_AM)),
                        otu_table(otu_table(AM_final, taxa_are_rows=FALSE)),
                        sample_data(ps_AM))


# Getting to know your phyloseq data ####

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


####### 
# Note: There are 16 trees with <10 ASV's -> for due diligence should run the analyses without 
# them and see if it makes a difference. For now will retain. 
######


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


#############
# Save for untransformed final phyloseq 

saveRDS(ps_AM_final, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_func_subset_final_2025.RDS")

# This should be used for the differential abundance analyses later on 


#############
# Save 'matrix' which now has functional info matched up to OTU's and taxa
write.csv(matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_OTU_funcs_AM.csv")

#############


####### more data exploration #########

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


## Using robust clr transformation, which only transforms the non-zero values in the dataset, and is thus 
# less affected by high zero counts than standard clr is 

clr_AM <- decostand(asv, method = "rclr")

clr_AM <- as.data.frame(clr_AM)

# load the clr data back into a phyloseq object 
ps_AM_clr <- phyloseq(tax_table(tax_table(ps_AM_final)),
                     otu_table(otu_table(clr_AM, taxa_are_rows=FALSE)),
                     sample_data(sample_metadata))


#####################
# save phyloseq object with transformed data 
saveRDS(ps_AM_clr, file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")

## This is the phyloseq object to be used in all downstream analyses 
#####################

#####################Beta-diversity##########################################

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
pca_clr_AM = prcomp(clr_sites_AM[7:369], center = T, scale = F)

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
#groups dispersions are different (not homogenous) p = 0.004

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
#                             diff                    lwr                     upr                    p adj
# Northern-Andrews   0.18681506357030031040 -2.3560911185586763672  2.72972124569927698801 0.9975066796109566258366
# Southern-Andrews   0.48539949048187480685 -2.2652216870396046922  3.23602066800335430585 0.9676583240917973061102
# WFDP-Andrews      -3.30092044917789273484 -5.9901817627034859015 -0.61165913565229956816 0.0094283591509782826989
# Southern-Northern  0.29858442691157449644 -2.5731028094355701263  3.17027166325871911923 0.9930299741736001717385
# WFDP-Northern     -3.48773551274819304524 -6.3007049639855949863 -0.67476606151079154827 0.0085306434493431781974
# WFDP-Southern     -3.78631993965976754168 -6.7883762590373253332 -0.78426362028220975020 0.0071682941336593808401

# differences between WFDP-Andrews, WFDP-Northern, and WFDP-Southern 


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

# vegan::adonis2(formula = aitchison_AM ~ Site, data = sample_data, permutations = 999, method = "euclidean")
#           Df SumOfSqs      R2      F Pr(>F)    
# Model      3    600.1 0.05313 2.3753  0.001 ***
# Residual 127  10695.4 0.94687                  
# Total    130  11295.5 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# adonis result is significant but could be impacted by the heterogeneity in group variances observed 
# in the betadispersion 


permanova.both <- vegan::adonis2(aitchison_AM ~ Site * Host_ID, data = sample_metadata, method = "euclidean", permutations = 999)
permanova.both

# explained 12.1%


permanova.both2 <- vegan::adonis2(aitchison_AM ~ Site + Host_ID, data = sample_metadata, method = "euclidean", permutations = 999)
permanova.both2

# explained 7.4%


##################

# Take distance to centroid values and regress with key environmental variables to assess if there are 
# relationships with environmental variation across the sites

# need to make dataset that has distance to centroid and the environmental variables 

# Have distance to centroid and environmental data for each individual tree 
distances_AM <- tibble::rownames_to_column(distances_AM, "Tree")

sample_metadata <- rownames_to_column(sample_metadata, "Tree") 

# pick some specific environmental variables of interest to compare 
AM_env2 <- dplyr::select(sample_metadata, Tree, elev, mean_precip_mm, mean_summer_precip_mm, MAT, pct_N, pct_C, Sand, 
                                          ph, org_matter, EC, avg_July_SPEI, count_mod_dry, count_sev_dry, apr1_SWE)

# combine with the distance to centroid summary table
centroid_enviro <- merge(distances_AM, AM_env2, by = "Tree")


### Testing these using the individual tree data

## Elevation
elev <- ggplot(centroid_enviro, aes(x = elev, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
elev

# test relationships 
lm_elev <- lm(DistanceToCentroid ~ elev, data = centroid_enviro)
summary(lm_elev) # not significant 


## Mean annual precip
MAP <- ggplot(centroid_enviro, aes(x = mean_precip_mm, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MAP

# test relationships 
lm_MAP <- lm(DistanceToCentroid ~ mean_precip_mm, data = centroid_enviro)
summary(lm_MAP) # SIGNIFICANT 

# Adjusted R-squared:  0.029659531203348743 p-value: 0.027950397405689189


## Mean summer precip
MSP <- ggplot(centroid_enviro, aes(x = mean_summer_precip_mm, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MSP

# test relationships 
lm_MSP <- lm(DistanceToCentroid ~ mean_summer_precip_mm, data = centroid_enviro)
summary(lm_MSP) # not significant 



## MAT
MAT <- ggplot(centroid_enviro, aes(x = MAT, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MAT

# test relationships 
lm_MAT <- lm(DistanceToCentroid ~ MAT, data = centroid_enviro)
summary(lm_MAT) # not significant 


## Soil percent C
pctC <- ggplot(centroid_enviro, aes(x = pct_C, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
pctC

# test relationships 
lm_pctC <- lm(DistanceToCentroid ~ pct_C, data = centroid_enviro)
summary(lm_pctC) # SIGNIFICANT 

# Adjusted R-squared:  0.04799103774521396 , p-value: 0.0070390795428453183


## Soil sand fraction
sand <- ggplot(centroid_enviro, aes(x = Sand, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
sand

# test relationships 
lm_sand <- lm(DistanceToCentroid ~ Sand, data = centroid_enviro)
summary(lm_sand) # SIGNIFICANT 

# Adjusted R-squared:  0.036548945655343834 , p-value: 0.016586243849719993


## pH
ph <- ggplot(centroid_enviro, aes(x = ph, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ph

# test relationships 
lm_ph <- lm(DistanceToCentroid ~ ph, data = centroid_enviro)
summary(lm_ph) # not significant 


## Cation exchange capacity (EC)
EC <- ggplot(centroid_enviro, aes(x = EC, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
EC

# test relationships 
lm_EC <- lm(DistanceToCentroid ~ EC, data = centroid_enviro)
summary(lm_EC) # SIGNIFICANT 

# Adjusted R-squared:  0.070076295061323068 ,  p-value: 0.001362559113534833


## Average July SPEI (avg_July_SPEI)
SPEI <- ggplot(centroid_enviro, aes(x = avg_July_SPEI, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SPEI

# test relationships 
lm_SPEI <- lm(DistanceToCentroid ~ avg_July_SPEI, data = centroid_enviro)
summary(lm_SPEI) # not significant 


## Count of Moderate Dry Months (count_mod_dry)
MOD <- ggplot(centroid_enviro, aes(x = count_mod_dry, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
MOD

# test relationships 
lm_MOD <- lm(DistanceToCentroid ~ count_mod_dry, data = centroid_enviro)
summary(lm_MOD) # not significant 

## Count of Severe Dry Months (count_mod_dry)
SEV <- ggplot(centroid_enviro, aes(x = count_sev_dry, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SEV

# test relationships 
lm_SEV <- lm(DistanceToCentroid ~ count_sev_dry, data = centroid_enviro)
summary(lm_SEV) # not significant 


## Average April 1 Snow Water Equivalent (apr1_SWE)
SWE <- ggplot(centroid_enviro, aes(x = apr1_SWE, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
SWE

# test relationships 
lm_SWE <- lm(DistanceToCentroid ~ apr1_SWE, data = centroid_enviro)
summary(lm_SWE) # Not significant 


### PLOT SIGNIFICANT RELATIONSHIPS ###

# Factors related to soil and drought 

#set colors for sites 
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# removed the legends from these so they will plot better as panels. Can add back in using theme(legend.position = "right")


# Mean Annual Precip
MAP_plot <- ggplot(centroid_enviro, aes(x = mean_precip_mm, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Mean Annual Precipitation (mm)", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

MAP_plot

# Soil C
soil_C_plot <- ggplot(centroid_enviro, aes(x = pct_C, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Soil C (%)", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                    name="Site",
                    breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                    labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

soil_C_plot

# Soil Sand Fraction
soil_Sand_plot <- ggplot(centroid_enviro, aes(x = Sand, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Soil Sand (%)", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

soil_Sand_plot


# Soil Cation Exchange Capacity
soil_EC_plot <- ggplot(centroid_enviro, aes(x = EC, y = DistanceToCentroid, colour = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "gray40") +
  theme_minimal() +
  labs(x = "Soil Cation Exchange Capacity (mmhos/cm)", y = "Distance to Centroid") +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  theme(legend.position = "NONE")

soil_EC_plot


# Could do a nice multi-panel figure with these using cowplot
cowplot::plot_grid(MAP_plot, soil_C_plot, soil_Sand_plot, soil_EC_plot, nrow = 2, ncol = 2)

##############

#check betadispersion between host tree taxa 
betadisper.host <- betadisper(aitchison_AM, sample_metadata$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper.host

permutest(betadisper.host)
#groups are different p = 0.006

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
#                    diff                     lwr                    upr                    p adj
# TABR-ALRU -3.31530020863627861161 -5.59337401409472967373 -1.0372264031778279936 0.0021750063351846371518
# THPL-ALRU -0.92395871263581419441 -3.14136766431316827308  1.2934502390415398843 0.5856714944069845163455
# THPL-TABR  2.39134149600046441719  0.18807643849146904458  4.5946065535094593457 0.0299318872859100082451

# TABR - ALRU 0.002
# TABR - THPL 0.030

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

# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# vegan::adonis2(formula = aitchison_AM ~ Host_ID, data = sample_metadata, permutations = 999, method = "euclidean")
#            Df              SumOfSqs                    R2                   F Pr(>F)
# Model      2   195.583621138572653 0.0185561788495433692 1.20059999999999989  0.131
# Residual 127 10344.497002377596800 0.9814438211504564746                           
# Total    129 10540.080623516170817 1.0000000000000000000  

# Not significant, explains only 1.9% 


########

# INTERPRETATION 

# There are significant differences in the beta dispersion of communities between sites and host tree taxa 
# The distance to centroid values for the tree communities are related to some environmental variables, but 
# simple linear relationships are an over-simplification here. 



# Save centroid_enviro_EM dataframe to use in future analyses 

write.csv(centroid_enviro, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_centroid_distance_enviro.csv")


