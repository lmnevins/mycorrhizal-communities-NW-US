# -----------------------------------------------------------------------------#
# Chapter 2: Merging AM fungal taxonomic data to functions and explore 
# Original Author: L. McKinley Nevins 
# July 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.2.1
#                     fundiversity v 1.1.1
#                     vegan v 2.6.10
#                     cluster v 2.1.8
#                     FD v 1.0.12.3
#                     phyloseq v 1.48.0
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(fundiversity); packageVersion("fundiversity")
library(vegan); packageVersion("vegan")
library(cluster); packageVersion("cluster")
library(FD); packageVersion("FD")
library(phyloseq); packageVersion("phyloseq")
#################################################################################
#                               Main workflow                                   #
#  Load in AM fungal taxa table, and the table of functions for genera that     #
#  had successful assignments. Match them up and format into matrices for       #
#  analyses of functional variation across sites and host species               #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/"
setwd(wd)

## Preparing the species x trait and species x site matrices 
###############

# Load in table of all Genera and OTUs matched with functional classifications 
funcs <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_OTU_funcs_AM.csv")


# Pull in phyloseq object with the clr transformed values for all downstream steps 
ps_AM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")


# reformat a tiny bit just to get OTU's as the species in rownames, and my traits only 
traits_AM <- funcs %>% dplyr::select(OTU2, Ancestral, Edaphophilic, Rhizophilic) %>% column_to_rownames(var = "OTU2") 

#OTUs are now the row names and there are columns for the binary coding of each trait level 

# Look at the reads of each tree
# Raw clr otu table 
otu_table <- otu_table(ps_AM_clr) %>% as.data.frame()


################## species by site matrix ########

## Keep the host tree code, because the host tree has it's own community and 
# the functioning likely varies between individual trees. 

# get clr abundance with each host tree as its own 'site' 

#######

trees_AM_clr <- otu_table(ps_AM_clr) %>% as("matrix") %>% as.data.frame()

#this already has each individual tree as it's own row 

# calling this "trees" because each tree is its own site 

####################
## 2. DATA ANALYSIS
# Relative Abundance of Functions 
####################

# I want to look at the traits for each tree community and compare the relative portions of traits that 
# are present using the clr values. 


# check structure 
str(trees_AM_clr)
str(traits_AM)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(trees_AM_clr)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_AM <- as.data.frame(trees_AM_clr)
tree_matrix_AM[] <- lapply(tree_matrix_AM, as.numeric)  
tree_matrix_AM <- as.matrix(tree_matrix_AM)
rownames(tree_matrix_AM) <- tree_ids 


# Keep rownames
otu_ids <- rownames(traits_AM)

traits_matrix_AM <- as.data.frame(traits_AM)
traits_matrix_AM[] <- lapply(traits_matrix_AM, as.numeric)
traits_matrix_AM <- as.matrix(traits_matrix_AM)
rownames(traits_matrix_AM) <- otu_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_AM)) # None        
sum(is.na(traits_matrix_AM))  # None   

# Make sure OTUs match between traits and trees
# Find shared OTUs between the two datasets 
shared_otus <- intersect(colnames(tree_matrix_AM), rownames(traits_matrix_AM))
# all are shared

# set the OTU's to be in the same order 
otu_order <- colnames(tree_matrix_AM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_AM <- traits_matrix_AM[otu_order, ]

# double check alignment 
all(colnames(tree_matrix_AM) == rownames(traits_matrix_AM)) #TRUE

# All match, good to go here 


check1 <- colSums(traits_matrix_AM) %>% as.data.frame() # All three guilds are present across samples
check2 <- rowSums(traits_matrix_AM) %>% as.data.frame() # each OTU is only assigned to one guild 




## Checking tree coverage 
# 1) Overlap check
common_otus <- intersect(colnames(tree_matrix_AM), rownames(traits_matrix_AM))
length(common_otus) # All OTUs have an assigned guild

# 2) Per-tree total *assigned* counts BEFORE CLR (use raw counts table)

# Read in raw final phyloseq 
ps_AM_raw <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_func_subset_final_2025.RDS")

raw_counts <- otu_table(ps_AM_raw) %>% as.data.frame()

raw_counts <- t(raw_counts)

assigned_counts_per_tree <- colSums(raw_counts[common_otus, , drop = FALSE]) %>% as.data.frame()

# 3) Per-tree assigned signal AFTER CLR
assigned_abs_sum <- rowSums(abs(tree_matrix_AM[, common_otus, drop = FALSE])) %>% as.data.frame()

# 4) How many assigned OTUs per tree?
assigned_otu_presence <- colSums(raw_counts[common_otus, , drop = FALSE] > 0)
summary(assigned_otu_presence)






############# 

# CLR-weighted trait composition: trees × traits
trait_sums_per_tree_AM <- tree_matrix_AM %*% traits_matrix_AM
# Output is a matrix of CLR-weighted totals of how much each trait is represented in a particular
# tree’s fungal community.


# More verbose way of doing this step that can be checked: 

# Make sure OTUs are aligned
tree_matrix_AM <- tree_matrix_AM[, rownames(traits_matrix_AM)]

# Initialize an empty results matrix
result <- matrix(0, nrow = nrow(tree_matrix_AM), ncol = ncol(traits_matrix_AM),
                 dimnames = list(rownames(tree_matrix_AM), colnames(traits_matrix_AM)))

# Loop over guilds
for (g in colnames(traits_matrix_AM)) {
  otus_in_g <- rownames(traits_matrix_AM)[traits_matrix_AM[, g] == 1]
  if (length(otus_in_g) > 0) {
    result[, g] <- rowSums(tree_matrix_AM[, otus_in_g, drop = FALSE])
  }
}


# Check which trees are losing all functional signal 
rowsums_per_tree <- rowSums(result) %>% as.matrix()
problem_trees <- names(rowsums_per_tree[rowsums_per_tree == 0]) %>% as.data.frame()



# Get dataframe for plotting 
trait_sums_df_AM <- as.data.frame(trait_sums_per_tree_AM)
trait_sums_df_AM$Sample_code <- rownames(trait_sums_df_AM)


# Reshape to long format
trait_sums_long_AM <- trait_sums_df_AM %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Sum")


# Join table of tree environmental data 
sample_data_AM <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_subset_sample_data_2025.csv")

sample_data_AM$Sample_code <- sample_data_AM$X

trait_sums_long_AM <- trait_sums_long_AM %>%
  left_join(sample_data_AM, by = "Sample_code")


#################

# Can now explore the clr-weighted trait profiles of each host tree! 

# set guild palette
              # Ancestral  Edaphophilic  Rhizophilic 
guild_colors <- c("darkgray", "brown4", "green4")


## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(trait_sums_long_AM$Sample_ID)
all_traits <- unique(trait_sums_long_AM$Trait)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
trait_sums_long_AM_full <- trait_sums_long_AM %>%
  dplyr::select(Sample_ID, Trait, CLR_Sum) %>%
  complete(Sample_ID = all_trees, Trait = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 
# Grab site from sample_data
site <- dplyr::select(sample_data_AM, Sample_ID, Site) 

trait_sums_long_AM_full <- merge(trait_sums_long_AM_full, site, by  = "Sample_ID")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

trait_sums_long_AM_full$Site <- factor(trait_sums_long_AM_full$Site, levels = site_order)

# and hosts
hosts <- dplyr::select(sample_data_AM, Sample_ID, Host_ID)

# merge hosts to data 
trait_sums_long_AM_full <- merge(trait_sums_long_AM_full, hosts, by = "Sample_ID")



## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos_AM <- trait_sums_long_AM_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(CLR_Sum)) %>%
  ungroup()

traits_neg_AM <- trait_sums_long_AM_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup() %>%
  mutate(CLR_Sum_scaled = -abs(CLR_Sum_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
traits_diverging_AM <- bind_rows(traits_pos_AM, traits_neg_AM)


# Diverging bar plot for trait representation in individual trees 

guild_diverging <- ggplot(traits_diverging_AM, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x") +
  theme_minimal() +
  labs(x = "Host Species",
       y = "Relative Guild Abundance",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) +
  geom_hline(yintercept = 0, color = "black")

guild_diverging

# Save traits_diverging_AM dataframe to use in future analyses 

write.csv(traits_diverging_AM, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_guild_clr_abund.csv")



#### 
## Calculate averages for each host in each site 

trait_avg_AM <- trait_sums_long_AM_full %>%
  group_by(Host_ID, Site, Trait) %>%
  summarise(mean_clr = mean(CLR_Sum, na.rm = TRUE), .groups = "drop")

#Split and scale positive and negative values by host
pos_avg_AM <- trait_avg_AM %>%
  filter(mean_clr > 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(mean_clr)) %>%
  ungroup()

neg_avg_AM <- trait_avg_AM %>%
  filter(mean_clr < 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(abs(mean_clr))) %>%
  ungroup() %>%
  mutate(mean_clr_scaled = -abs(mean_clr_scaled))  # ensure values are negative

trait_avg_scaled_AM <- bind_rows(pos_avg_AM, neg_avg_AM)

# Save traits_avg_scaled dataframe to use in future analyses 
write.csv(trait_avg_scaled_AM, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_guild_avg_scaled.csv")



# PLOT

host_guild_plot_AM <- ggplot(trait_avg_scaled_AM, aes(x = Host_ID, y = mean_clr_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme_minimal() +
  labs(y = "Relative Guild Abundance",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  theme(
    axis.text.x = element_text(size = 10, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 12)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"))

host_guild_plot_AM



########################################################

# Get weighted taxon abundance in tree communities 

# Get dataframe for plotting 
tree_abund_df <- as.data.frame(tree_matrix_AM) %>% rownames_to_column(var = "Sample_code")


# Join table of tree environmental data 

#subset sample data 
sample_data <- dplyr::select(sample_data_AM, Sample_code, Site, Host_ID)


tree_abund_full <- tree_abund_df %>%
  merge(sample_data, by = "Sample_code")


# Reshape to long format
tree_abund_long <- tree_abund_full %>%
  pivot_longer(-c(Sample_code, Site, Host_ID), names_to = "OTU", values_to = "CLR_Abund")


## Merge with taxa data 
# Pull out tax table from the trimmed table for taxa that have functional assignments 
AM_tax <- tax_table(ps_AM_clr) %>% as("matrix") %>% as.data.frame()

AM_tax <- as.data.frame(AM_tax) %>% rownames_to_column(var = "OTU")

tree_tax_full <- merge(AM_tax, tree_abund_long, by = "OTU")


# Clean up tree names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Tree = Sample_code %>%
           str_replace_all("_FWD_filt.fastq.gz", ""))


#### manipulations to plot 

# plot per tree just to get a look at things 
taxa_tree_plot <- ggplot(tree_tax_full, aes(x = Tree, y = CLR_Abund, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "CLR-Weighted Fungal Taxonomic Composition per Tree",
    x = "",
    y = "Relative Family Abundance",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 9)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

taxa_tree_plot



## Fill in some trees that have no bars 
# Get list of all trees
all_trees <- unique(tree_tax_full$Sample_code)


# Expand to all Trees and fill missing with 0
# Putting a pin in this, not sure I need to do this 
tree_tax_plotting <- tree_tax_full %>%
  dplyr::select(Sample_code, Tree, CLR_Abund, Genus, Family) %>%
  complete(Sample_code = all_trees, fill = list(CLR_Abund = 0))

# Add a couple other data columns for plotting later 
# Grab site from sample_data
site <- dplyr::select(sample_data, Sample_code, Site) 

tree_tax_plotting <- merge(tree_tax_plotting, site, by  = "Sample_code")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

tree_tax_plotting$Site <- factor(tree_tax_plotting$Site, levels = site_order)

# and hosts
hosts <- dplyr::select(sample_data, Sample_code, Host_ID)

# merge hosts to data 
tree_tax_plotting <- merge(tree_tax_plotting, hosts, by = "Sample_code")


## Here these can diverge to consider the taxa that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
tree_taxa_pos <- tree_tax_plotting %>%
  filter(CLR_Abund > 0) %>%
  group_by(Sample_code) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(CLR_Abund)) %>%
  ungroup()

tree_taxa_neg <- tree_tax_plotting %>%
  filter(CLR_Abund < 0) %>%
  group_by(Sample_code) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup() %>%
  mutate(CLR_Abund_scaled = -abs(CLR_Abund_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
tree_taxa_diverging <- bind_rows(tree_taxa_pos, tree_taxa_neg)


# Diverging bar plot for taxon representation in individual trees 

tree_taxa_diverging <- ggplot(tree_taxa_diverging, aes(x = Tree, y = CLR_Abund_scaled, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Tree",
    y = "Relative Genus Abundance",
    fill = "Genus") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 9)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12)) +
  geom_hline(yintercept = 0, color = "black")

tree_taxa_diverging



## Calculate averages for each host in each site 

tree_taxa_avg <- tree_tax_plotting %>%
  group_by(Host_ID, Site, Genus, Family) %>%
  summarise(mean_clr = mean(CLR_Abund, na.rm = TRUE), .groups = "drop")

#Split and scale positive and negative values by host
tree_taxa_pos_avg <- tree_taxa_avg %>%
  filter(mean_clr > 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(mean_clr)) %>%
  ungroup()

tree_taxa_neg_avg <- tree_taxa_avg %>%
  filter(mean_clr < 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(abs(mean_clr))) %>%
  ungroup() %>%
  mutate(mean_clr_scaled = -abs(mean_clr_scaled))  # ensure values are negative

tree_taxa_avg_scaled <- bind_rows(tree_taxa_pos_avg, tree_taxa_neg_avg)


# PLOT

host_tree_taxa_plot <- ggplot(tree_taxa_avg_scaled, aes(x = Host_ID, y = mean_clr_scaled, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  theme_minimal() +
  labs(
    x = "Tree",
    y = "Relative Genus Abundance",
    fill = "Genus") +
  theme_minimal(base_size = 12) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 10)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12)) +
  geom_hline(yintercept = 0, color = "black")

host_tree_taxa_plot


# average for families instead 

tree_taxa_avg_fam <- tree_tax_plotting %>%
  group_by(Host_ID, Site, Family) %>%
  summarise(mean_clr = mean(CLR_Abund, na.rm = TRUE), .groups = "drop")

#Split and scale positive and negative values by host
tree_taxa_pos_avg_fam <- tree_taxa_avg_fam %>%
  filter(mean_clr > 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(mean_clr)) %>%
  ungroup()

tree_taxa_neg_avg_fam <- tree_taxa_avg_fam %>%
  filter(mean_clr < 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(abs(mean_clr))) %>%
  ungroup() %>%
  mutate(mean_clr_scaled = -abs(mean_clr_scaled))  # ensure values are negative

tree_taxa_avg_scaled_fam <- bind_rows(tree_taxa_pos_avg_fam, tree_taxa_neg_avg_fam)


# PLOT

host_tree_taxa_fam_plot <- ggplot(tree_taxa_avg_scaled_fam, aes(x = Host_ID, y = mean_clr_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  theme_minimal() +
  labs(,
    y = "Relative Family Abundance",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 12, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1)

host_tree_taxa_fam_plot


########################################################################################

#####################
## 3. DATA ANALYSIS
# PERMANOVA for traits x host 
# Beta dispersion of functions
####################

# Grab the scaled average CLR-weighted abundance values of each guild for each host species across sites 
mod_traits <- trait_avg_scaled_AM

# Create a site-Host ID column for merging 
mod_traits <- mod_traits %>% mutate(id_code = paste(Site, Host_ID, sep = "_"))

# Trim site column so it is not duplicated when merging later
mod_traits <- dplyr::select(mod_traits, Trait, mean_clr, mean_clr_scaled, id_code, Host_ID)


# Load in data for the environmental variables that were significant from forward selection
# which was ph and count_mod_dry for the AM communities 

# sample_data_AM has this info 
dat <- dplyr::select(sample_data_AM, Host_ID, Site, Sample_ID, ph, count_mod_dry)

# This is for individual trees, so take averages of the environmental variables to get them 
# at the species-site level 
dat <- dat %>%
  group_by(Site, Host_ID) %>%
  summarize(ph = mean(ph), count_mod_dry = mean(count_mod_dry))


# Create a site-Host ID column for merging 
dat <- dat %>% mutate(id_code = paste(Site, Host_ID, sep = "_"))

# Trim Host_ID column so they are not duplicated when merging 
dat <- dplyr::select(dat, ph, count_mod_dry, id_code)


# Merge dataframes by id_code column 
mod_data <- merge(dat, mod_traits, by = "id_code")


## Dataset now has the environmental variables for each host at each site, 
# and the mean CLR scaled abundance for each of the guilds for each host species in each site 


# Also make a dataset for each individual tree instead 
indiv_traits <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_guild_clr_abund.csv")

# reformat to get abundance for each trait and the sample_ID as the rownames 

indiv_matrix <- indiv_traits %>%
  mutate(CLR_Sum_scaled = replace_na(CLR_Sum_scaled, 0)) %>%
  pivot_wider(
    id_cols = Sample_ID,
    names_from = Trait,
    values_from = CLR_Sum_scaled,
    values_fill = 0  
  ) %>%
  column_to_rownames("Sample_ID")

## now this dataset is a trait matrix of the CLR weighted abundance values for each guild present 
# in each tree community. Now I can analyze the variation at the individual tree level and how the functional 
# variation clusters 

# save indiv_matrix
write.csv(indiv_matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_tree_guild_clr.csv")


# Get metadata for the individual tree level 
indiv_meta <- indiv_traits %>% 
  dplyr::select(Site, Host_ID, Sample_ID) 


indiv_meta <- indiv_meta %>% distinct() %>% column_to_rownames("Sample_ID") 

indiv_meta_data <- rownames_to_column(indiv_meta, var = "Sample_ID")


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
pca_traits_AM = prcomp(traits_for_pca_AM[4:6], center = T, scale = F)

sd.pca_traits_AM = pca_traits_AM$sdev
loadings.pca_traits_AM = pca_traits_AM$rotation
names.pca_traits_AM = colnames(traits_for_pca_AM[4:6])
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
betadisper_traits_site <- vegan::betadisper(aitchison_AM_traits, indiv_meta_data$Site, type = "median", sqrt.dist = FALSE)

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
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      3  1.055 0.35172 0.9372    999  0.414
# Residuals 101 37.905 0.37530 



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
#                        diff        lwr       upr     p adj
# Northern-Andrews   0.18865710 -0.2060459 0.5833601 0.5975892
# Southern-Andrews  -0.04800003 -0.4834170 0.3874170 0.9916297
# WFDP-Andrews       0.16112870 -0.3119752 0.6342326 0.8102095
# Southern-Northern -0.23665714 -0.6827847 0.2094704 0.5111634
# WFDP-Northern     -0.02752841 -0.5105078 0.4554510 0.9988158
# WFDP-Southern      0.20912873 -0.3076562 0.7259136 0.7162261


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
betadisper_traits_host <- betadisper(aitchison_AM_traits, indiv_meta_data$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper_traits_host

permutest(betadisper_traits_host)
#groups are not different p = 0.26

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      2  1.389 0.69448 1.3943    999   0.26
# Residuals 102 50.806 0.49809             


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
#                 diff        lwr       upr     p adj
# TABR-ALRU -0.28304965 -0.7157487 0.1496494 0.2695539
# THPL-ALRU -0.04061114 -0.4238980 0.3426757 0.9656093
# THPL-TABR  0.24243851 -0.1679199 0.6527969 0.3420506


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

# not significant, and only 3.1% explained


permanova.host <- adonis2(indiv_matrix ~ Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.host

# not significant, and only 0.9% explained

permanova.both <- vegan::adonis2(indiv_matrix ~ Site * Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both

# not significant, and only 6.5% explained


permanova.both2 <- vegan::adonis2(indiv_matrix ~ Site + Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both2

# not significant, and only 3.8% explained


########################################################################################


