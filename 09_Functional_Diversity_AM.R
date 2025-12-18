# -----------------------------------------------------------------------------#
# Chapter 2: Merging AM fungal taxonomic data to functions and explore 
# Original Author: L. McKinley Nevins 
# December 16, 2025
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

# Load in table of all Genera and ASVs matched with functional classifications 
funcs <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_ASV_funcs_AM.csv")


# Pull in phyloseq object with the clr transformed values for all downstream steps 
ps_AM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_transformed_final.RDS")


# Pull in phyloseq object with the raw ASV abundances 
ps_AM_raw <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_func_subset_final_2025.RDS")


# reformat a tiny bit just to get ASVs as the species in rownames, and my traits only 
traits_AM <- funcs %>% dplyr::select(ASV2, Ancestral, Edaphophilic, Rhizophilic) %>% column_to_rownames(var = "ASV2") 

#ASVs are now the row names and there are columns for the binary coding of each trait level 

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


# get raw ASVs with each host tree as its own 'site' 

trees_AM_raw <- otu_table(ps_AM_raw) %>% as("matrix") %>% as.data.frame()

####################
## 2. DATA ANALYSIS
# Relative Abundance of Functions 
####################


## Want to use the raw ASV abundances first, then do clr transformation of the traits after they have been 
# aggregated across the ASVs

# This is an edited version of what is down below 


# I want to look at the traits for each tree community and compare the relative portions of traits that 
# are present using the clr values. 


# check structure 
str(trees_AM_raw)
str(traits_AM)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(trees_AM_raw)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_AM <- as.data.frame(trees_AM_raw)
tree_matrix_AM[] <- lapply(tree_matrix_AM, as.numeric)  
tree_matrix_AM <- as.matrix(tree_matrix_AM)
rownames(tree_matrix_AM) <- tree_ids 


# Keep rownames
asv_ids <- rownames(traits_AM)

traits_matrix_AM <- as.data.frame(traits_AM)
traits_matrix_AM[] <- lapply(traits_matrix_AM, as.numeric)
traits_matrix_AM <- as.matrix(traits_matrix_AM)
rownames(traits_matrix_AM) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_AM)) # None        
sum(is.na(traits_matrix_AM))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_AM), rownames(traits_matrix_AM))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_AM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_AM <- traits_matrix_AM[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_AM) == rownames(traits_matrix_AM)) #TRUE

# All match, good to go here 


check1 <- colSums(traits_matrix_AM) %>% as.data.frame() # All three guilds are present across samples
check2 <- rowSums(traits_matrix_AM) %>% as.data.frame() # each ASV is only assigned to one guild 


############# 

# Trait composition: trees × traits using raw ASV counts 

# First for exploration types 

trait_abund_per_tree <- tree_matrix_AM %*% traits_matrix_AM
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion, which makes it 
# compatible with aitchison distance that uses compositions 
trait_prop_per_tree <- trait_abund_per_tree / rowSums(trait_abund_per_tree)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
trait_clr_per_tree <- decostand(trait_prop_per_tree, method = "clr", pseudocount = 1e-06)

## Save file of CLR abundance for each trait 

write.csv(trait_clr_per_tree, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_trait_clr_per_tree.csv")


## After this point the variable names are the same as the other option down below 
# Get dataframe for plotting 
trait_sums_df_AM <- as.data.frame(trait_clr_per_tree)
trait_sums_df_AM$Sample_code <- rownames(trait_sums_df_AM)


# Reshape to long format
trait_sums_long_AM <- trait_sums_df_AM %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Sum")


# Join table of tree environmental data 
sample_data <- sample_data(ps_AM_raw) %>% as("matrix") %>% as.data.frame()

sample_data <- sample_data %>% rownames_to_column(var = "Sample_code")

trait_sums_long_AM <- trait_sums_long_AM %>%
  left_join(sample_data, by = "Sample_code")


### SKIP: ORIGINAL VERSION WITH CLR-WEIGHTED RELATIVE ABUNDANCE #### 

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
site <- dplyr::select(sample_data, Sample_ID, Site) 

trait_sums_long_AM_full <- merge(trait_sums_long_AM_full, site, by  = "Sample_ID")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

trait_sums_long_AM_full$Site <- factor(trait_sums_long_AM_full$Site, levels = site_order)

# and hosts
hosts <- dplyr::select(sample_data, Sample_ID, Host_ID)

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
  ungroup()

# Merge for plotting 

# Combine
traits_diverging_AM <- bind_rows(traits_pos_AM, traits_neg_AM)


# Diverging bar plot for trait representation in individual trees 

guild_diverging <- ggplot(traits_diverging_AM, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x",nrow = 4) +
  theme_minimal() +
  labs(x = "Host Species",
       y = "Relative Over- or Under-Representation of Guilds",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"))

guild_diverging

# Save traits_diverging_AM dataframe to use in future analyses 

write.csv(traits_diverging_AM, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_guild_clr_abund.csv")


# Subset to one location and facet by host species 

guild_north <- dplyr::filter(traits_diverging_AM, Site == "Northern")


guild_north_plot <- ggplot(guild_north, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(x = "",
       y = "Relative Over- or Under-Representation of Guilds",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

guild_north_plot


guild_wfdp <- dplyr::filter(traits_diverging_AM, Site == "WFDP")


guild_wfdp_plot <- ggplot(guild_wfdp, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(x = "",
       y = "Relative Over- or Under-Representation of Guilds",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

guild_wfdp_plot


guild_andrews <- dplyr::filter(traits_diverging_AM, Site == "Andrews")


guild_andrews_plot <- ggplot(guild_andrews, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(x = "",
       y = "Relative Over- or Under-Representation of Guilds",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

guild_andrews_plot


guild_south <- dplyr::filter(traits_diverging_AM, Site == "Southern")


guild_south_plot <- ggplot(guild_south, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(x = "",
       y = "Relative Over- or Under-Representation of Guilds",
       fill = "Guild") +
  scale_fill_manual(values=guild_colors, 
                    name="Guild",
                    breaks=c("Ancestral", "Edaphophilic", "Rhizophilic"),
                    labels=c("Ancestral", "Edaphophilic", "Rhizophilic")) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=11, face="bold"), 
        legend.position = "bottom")

guild_south_plot



########################################################

# Get weighted taxon abundance in tree communities 

## Using the clr transformed values from the trees_AM_clr object, pulled from the transformed phyloseq

# Keep rownames
tree_ids <- rownames(trees_AM_clr)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_AM_clr <- as.data.frame(trees_AM_clr)
tree_matrix_AM_clr[] <- lapply(tree_matrix_AM_clr, as.numeric)  
tree_matrix_AM_clr <- as.matrix(tree_matrix_AM_clr)
rownames(tree_matrix_AM_clr) <- tree_ids 


# Get dataframe for plotting 
tree_abund_df <- as.data.frame(tree_matrix_AM_clr) %>% rownames_to_column(var = "Sample_code")


# Join table of tree environmental data 

#subset sample data 
sample_data <- dplyr::select(sample_data, Sample_code, Site, Host_ID)


tree_abund_full <- tree_abund_df %>%
  merge(sample_data, by = "Sample_code")


# Reshape to long format
tree_abund_long <- tree_abund_full %>%
  pivot_longer(-c(Sample_code, Site, Host_ID), names_to = "ASV", values_to = "CLR_Abund")


## Merge with taxa data 
# Pull out tax table from the trimmed table for taxa that have functional assignments 
AM_tax <- tax_table(ps_AM_clr) %>% as("matrix") %>% as.data.frame()

AM_tax <- as.data.frame(AM_tax) %>% rownames_to_column(var = "ASV")

tree_tax_full <- merge(AM_tax, tree_abund_long, by = "ASV")


# Clean up tree names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Tree = Sample_code %>%
           str_replace_all("_FWD_filt.fastq.gz", ""))


#### manipulations to plot 

# plot per tree just to get a look at things - not scaled yet 
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
  ungroup()

# Merge for plotting 

# Combine
tree_taxa_diverging <- bind_rows(tree_taxa_pos, tree_taxa_neg)


# Diverging bar plot for taxon representation in individual trees 

tree_taxa_plot <- ggplot(tree_taxa_diverging, aes(x = Tree, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Tree",
    y = "Relative Over- or Under-Representation of Arbuscular Mycorrhizal Families",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 8)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 0, colour="black")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 11)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 11)) +
  theme(legend.title = element_text(colour="black", size=11, face = "bold"), 
        legend.position = "bottom")


tree_taxa_plot




# Subset to one location and facet by host species 

fam_north <- dplyr::filter(tree_taxa_diverging, Site == "Northern")

tree_fam_north_plot <- ggplot(fam_north, aes(x = Tree, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Arbuscular Mycorrhizal Families",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 8)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 0, colour="black")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 11)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 11)) +
  theme(legend.title = element_text(colour="black", size=11, face = "bold"), 
        legend.position = "none")

tree_fam_north_plot


fam_wfdp <- dplyr::filter(tree_taxa_diverging, Site == "WFDP")

tree_fam_wfdp_plot <- ggplot(fam_wfdp, aes(x = Tree, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Arbuscular Mycorrhizal Families",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 8)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 0, colour="black")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 11)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 11)) +
  theme(legend.title = element_text(colour="black", size=11, face = "bold"), 
        legend.position = "none")

tree_fam_wfdp_plot


fam_andrews <- dplyr::filter(tree_taxa_diverging, Site == "Andrews")

tree_fam_andrews_plot <- ggplot(fam_andrews, aes(x = Tree, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Arbuscular Mycorrhizal Families",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 8)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 0, colour="black")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 11)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 11)) +
  theme(legend.title = element_text(colour="black", size=11, face = "bold"), 
        legend.position = "none")

tree_fam_andrews_plot


fam_south <- dplyr::filter(tree_taxa_diverging, Site == "Southern")


tree_fam_south_plot <- ggplot(fam_south, aes(x = Tree, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Arbuscular Mycorrhizal Families",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 8)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 0, colour="black")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 11)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 11)) +
  theme(legend.title = element_text(colour="black", size=11, face = "bold"), 
        legend.position = "none")

tree_fam_south_plot



## Done visualizing the functional and taxonomic variation across the sites and hosts 

## -- END -- ## 

