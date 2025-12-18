# -----------------------------------------------------------------------------#
# Chapter 2: Explore and visualize EM community functional variation 
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
#                     ade4 v 1.7.23
#                     phyloseq v 1.48.0
#                     ape v 5.8.1
#                     ggh4x v 0.3.1.9000
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(fundiversity); packageVersion("fundiversity")
library(vegan); packageVersion("vegan")
library(cluster); packageVersion("cluster")
library(FD); packageVersion("FD")
library(ade4); packageVersion("ade4")
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
library(ggh4x); packageVersion("ggh4x")

#################################################################################
#                               Main workflow                                   #
#  Load in EM fungal taxa table, and the table of functions for genera that     #
#  had successful assignments. Match them up and format into matrices for       #
#  analyses of functional variation across sites and host species, including    #
#  7 EMF hosts: ABAM, ABGR, ABPR, ALRU, PSME, TABR, TSHE. Visualize variation   #
#  in exploration types, hydrophobic traits, and fungal families across sites   #
#  and host tree species. 
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
funcs <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_OTU_funcs_EM_no_THPL.csv")


# this has a categorical 'hydro' column for hydrophobicity, and a 'hydro_binary' column where it's been coded 
# as hydrophilic = 1 and hydrophobic = 0.

# The Exploration Types have all been one-hot encoded as separate binary columns - there were 14 categories


# Pull in phyloseq object with the clr transformed values 
ps_EM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final_no_THPL.RDS")

# 371 trees and 1,627 ASVs


# Pull in phyloseq object with the raw ASV abundances 
ps_EM_raw <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_func_final_no_THPL.RDS")


# reformat funcs a tiny bit just to get ASVs as the species in rownames, and my traits only 
traits_EM_hydro <- funcs %>% dplyr::select(ASV2, hydrophilic, hydrophobic) %>% column_to_rownames(var = "ASV2") 

traits_EM_ET <- funcs %>% dplyr::select(ASV2, ET_contact, ET_contact_short, ET_short, ET_contact_medium,
                                           ET_contact_medium_fringe, ET_contact_medium_smooth, ET_medium_smooth, 
                                           ET_medium_fringe, ET_medium_mat, ET_medium_long, ET_medium_long_smooth,
                                           ET_medium_long_fringe, ET_contact_long_smooth, ET_long) %>% column_to_rownames(var = "ASV2") 

#ASVs are now the row names and there are columns for the binary coding of each trait level, with separate datasets for 
# hydrophobic traits and exploration types 


################## species by site matrix ########

## Keep the host tree code, because the host tree has it's own community and 
# the functioning likely varies between individual trees. 

# get clr abundance with each host tree as its own 'site' 

#######

trees_EM_clr <- otu_table(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

#this already has each individual tree as it's own row 

# calling this "trees" because each tree is its own site 


# get raw ASVs with each host tree as its own 'site' 

trees_EM_raw <- otu_table(ps_EM_raw) %>% as("matrix") %>% as.data.frame()

####################
## 2. DATA ANALYSIS
# Relative Abundance of Functions 
####################

## Want to use the raw ASV abundances first, then do clr transformation of the traits after they have been 
# aggregated across the ASVs

# This is an edited version of what is down below 

### Need to do ET and hydro traits separately, they cannot be considered proportions
# of the same whole 

# Do for Exploration types first 

# check structure 
str(trees_EM_raw)
str(traits_EM_ET)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(trees_EM_raw)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_EM <- as.data.frame(trees_EM_raw)
tree_matrix_EM[] <- lapply(tree_matrix_EM, as.numeric)  
tree_matrix_EM <- as.matrix(tree_matrix_EM)
rownames(tree_matrix_EM) <- tree_ids 


# Keep rownames
asv_ids <- rownames(traits_EM_ET)

traits_matrix_EM_ET <- as.data.frame(traits_EM_ET)
traits_matrix_EM_ET[] <- lapply(traits_matrix_EM_ET, as.numeric)
traits_matrix_EM_ET <- as.matrix(traits_matrix_EM_ET)
rownames(traits_matrix_EM_ET) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_EM)) # None        
sum(is.na(traits_matrix_EM_ET))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM_ET))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM_ET <- traits_matrix_EM_ET[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM_ET)) #TRUE

# All match, good to go here 


# Repeat steps for the hydrophobic traits matrix 

str(traits_EM_hydro)

# make numeric matrix

# Keep rownames
asv_ids <- rownames(traits_EM_hydro)

traits_matrix_EM_hydro <- as.data.frame(traits_EM_hydro)
traits_matrix_EM_hydro[] <- lapply(traits_matrix_EM_hydro, as.numeric)
traits_matrix_EM_hydro <- as.matrix(traits_matrix_EM_hydro)
rownames(traits_matrix_EM_hydro) <- asv_ids

# check for any NAs in the data 
sum(is.na(traits_matrix_EM_hydro))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM_hydro))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM_hydro <- traits_matrix_EM_hydro[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM_hydro)) #TRUE

# All match, good to go here 


############# 

# Trait composition: trees × traits using raw ASV counts 

# First for exploration types 

trait_abund_per_tree_ET <- tree_matrix_EM %*% traits_matrix_EM_ET
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion, which makes it 
# compatible with aitchison distance that uses compositions 
trait_prop_per_tree_ET <- trait_abund_per_tree_ET / rowSums(trait_abund_per_tree_ET)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
# Using the same clr and the same pseudocount as the taxonomic analyses 
trait_clr_per_tree_ET <- decostand(trait_prop_per_tree_ET, method = "clr", pseudocount = 1e-06)


## Save file of CLR abundance for each trait 

write.csv(trait_clr_per_tree_ET, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_trait_clr_per_tree.csv")


## After this point the variable names are the same as the other option down below 
# Get dataframe for plotting 
trait_sums_df_ET <- as.data.frame(trait_clr_per_tree_ET)
trait_sums_df_ET$Sample_code <- rownames(trait_sums_df_ET)


# Reshape to long format
trait_sums_long_ET <- trait_sums_df_ET %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Sum")


# Join table of tree environmental data 
sample_data <- sample_data(ps_EM_raw) %>% as("matrix") %>% as.data.frame()

sample_data <- sample_data %>% rownames_to_column(var = "Sample_code")

trait_sums_long_ET <- trait_sums_long_ET %>%
  left_join(sample_data, by = "Sample_code")


# And repeat steps for the hydrophobic traits 

trait_abund_per_tree_hydro <- tree_matrix_EM %*% traits_matrix_EM_hydro
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion, which makes it 
# compatible with aitchison distance that uses compositions 
trait_prop_per_tree_hydro <- trait_abund_per_tree_hydro / rowSums(trait_abund_per_tree_hydro)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
trait_clr_per_tree_hydro <- decostand(trait_prop_per_tree_hydro, method = "clr", pseudocount = 1e-06)


## After this point the variable names are the same as the other option down below 
# Get dataframe for plotting 
trait_sums_df_hydro <- as.data.frame(trait_clr_per_tree_hydro)
trait_sums_df_hydro$Sample_code <- rownames(trait_sums_df_hydro)


# Reshape to long format
trait_sums_long_hydro <- trait_sums_df_hydro %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Sum")


# Join table of tree environmental data 
trait_sums_long_hydro <- trait_sums_long_hydro %>%
  left_join(sample_data, by = "Sample_code")


### SKIP: ORIGINAL VERSION WITH CLR-WEIGHTED RELATIVE ABUNDANCE #### 


# I want to look at the traits for each tree community and compare the relative portions of traits that 
# are present using the clr values. 

# check structure 
str(trees_EM_clr)
str(traits_EM)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(trees_EM_clr)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_EM <- as.data.frame(trees_EM_clr)
tree_matrix_EM[] <- lapply(tree_matrix_EM, as.numeric)  
tree_matrix_EM <- as.matrix(tree_matrix_EM)
rownames(tree_matrix_EM) <- tree_ids 


# Keep rownames
asv_ids <- rownames(traits_EM)

traits_matrix_EM <- as.data.frame(traits_EM)
traits_matrix_EM[] <- lapply(traits_matrix_EM, as.numeric)
traits_matrix_EM <- as.matrix(traits_matrix_EM)
rownames(traits_matrix_EM) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_EM)) # None        
sum(is.na(traits_matrix_EM))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM <- traits_matrix_EM[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM)) #TRUE

# All match, good to go here 

############ -- 

# CLR-weighted trait composition: trees × traits
trait_sums_per_tree <- tree_matrix_EM %*% traits_matrix_EM
# Output is a matrix of CLR-weighted totals of how much each trait is represented in a particular
# tree’s fungal community.


# Get dataframe for plotting 
trait_sums_df <- as.data.frame(trait_sums_per_tree)
trait_sums_df$Sample_code <- rownames(trait_sums_df)


# Reshape to long format
trait_sums_long <- trait_sums_df %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Sum")


# Join table of tree environmental data 
sample_data <- sample_data(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

sample_data <- sample_data %>% rownames_to_column(var = "Sample_code")

trait_sums_long <- trait_sums_long %>%
  left_join(sample_data, by = "Sample_code")


#################


# This is now using the new method for calculating the clr abundance of each trait in each 
# host tree community 


# Can now explore the clr-weighted trait profiles of each host tree! 

# Reformat names in the exploraiton type and hydrophobicity datasets to make them a bit nicer 

# Clean up the trait names 
trait_sums_long_ET <- trait_sums_long_ET %>%
  mutate(Trait_clean = Trait %>%
           str_remove("^ET_") %>%                # remove prefix
           str_replace_all("_", "-") %>%         # underscores → dashes
           str_to_title()                        # capitalize each word
  )


# Hyphal hydrophobicity
trait_sums_long_hydro <- trait_sums_long_hydro %>%
  mutate(Trait_clean = Trait %>%
           str_to_title())



# Set exploration type order to reflect range of short to long distance investment 
exploration_order <- c("Contact", "Contact-Short", "Short", "Contact-Medium", "Contact-Medium-Smooth",
                       "Contact-Medium-Fringe", "Contact-Long-Smooth", "Medium-Smooth", "Medium-Fringe", 
                       "Medium-Mat", "Medium-Long", "Medium-Long-Smooth", "Medium-Long-Fringe", "Long")


# apply order to the trait column 
trait_sums_long_ET <- trait_sums_long_ET %>%
  mutate(Trait_clean = factor(Trait_clean, levels = exploration_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(14, option = "D", direction = -1)
print(viridis_colors)


# plot per tree 
ET_tree_plot <- ggplot(trait_sums_long_ET, aes(x = Sample_ID, y = CLR_Sum, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative Over- and Under-Representation of Exploration Types",
       fill = "Exploration Type"
  ) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

ET_tree_plot


## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(trait_sums_long_ET$Sample_ID)
all_traits <- unique(trait_sums_long_ET$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
trait_sums_long_ET_full <- trait_sums_long_ET %>%
  dplyr::select(Sample_ID, Trait_clean, CLR_Sum) %>%
  complete(Sample_ID = all_trees, Trait_clean = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 
# Grab site from sample_data
site <- dplyr::select(sample_data, Sample_ID, Site) 

trait_sums_long_ET_full <- merge(trait_sums_long_ET_full, site, by  = "Sample_ID")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

trait_sums_long_ET_full$Site <- factor(trait_sums_long_ET_full$Site, levels = site_order)

# and hosts
hosts <- dplyr::select(sample_data, Sample_ID, Host_ID)

# merge hosts to data 
trait_sums_long_ET_full <- merge(trait_sums_long_ET_full, hosts, by = "Sample_ID")


# Clean sample_ID for less visual clutter 
trait_sums_long_ET_full <- trait_sums_long_ET_full %>%
  mutate(Tree_ID = Sample_ID %>%
           str_remove("^A-") %>%
           str_remove("^N-") %>%
           str_remove("^S-") %>%
           str_remove("^W-"))

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos <- trait_sums_long_ET_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup()

traits_neg <- trait_sums_long_ET_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup()

# Merge for plotting 

# Combine
traits_diverging <- bind_rows(traits_pos, traits_neg)


# Make a couple variables factors 
traits_diverging$Host_ID <- as.factor(traits_diverging$Host_ID)
traits_diverging$Tree_ID <- as.factor(traits_diverging$Tree_ID)



# Diverging bar plot for trait representation in individual trees 

ET_diverging <- ggplot(traits_diverging, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x", nrow = 4) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal() +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Exploration Types",
    fill = "Exploration Type") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 7, colour="black"),
  axis.text.y = element_text(size = 11, colour="black"),
  axis.title.y = element_text(size = 12, colour="black"),
  legend.text = element_text(size = 11, colour="black"),
  strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"))

ET_diverging



# Subset to one location and facet by host species 

ET_north <- dplyr::filter(traits_diverging, Site == "Northern")


ET_north_plot <- ggplot(ET_north, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Exploration Types",
    fill = "Exploration Type") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

ET_north_plot


ET_wfdp <- dplyr::filter(traits_diverging, Site == "WFDP")


ET_wfdp_plot <- ggplot(ET_wfdp, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Exploration Types",
    fill = "Exploration Type") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

ET_wfdp_plot


ET_andrews <- dplyr::filter(traits_diverging, Site == "Andrews")


ET_andrews_plot <- ggplot(ET_andrews, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Exploration Types",
    fill = "Exploration Type") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

ET_andrews_plot


ET_south <- dplyr::filter(traits_diverging, Site == "Southern")


ET_south_plot <- ggplot(ET_south, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Exploration Types",
    fill = "Exploration Type") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 10, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=11, face="bold"), 
        legend.position = "none")

ET_south_plot


#### Traits diverging dataframe has each individual tree, the site and host, and the CLR raw
# and CLR scaled weighted abundance for each level of the exploration type trait 

# Save traits_diverging dataframe to use in future analyses 

# This is saving the new way of calculating the clr-relative trait values 

write.csv(traits_diverging, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_clr_abund_no_THPL.csv")



# This shows there is a lot of variation within a given host tree species, so I want to demonstrate 
# this by keeping the individual trees separate. 


######## hydrophobicity analyses ######

# set hydro palette
# hydrophilic   hydrophobic
hydro_colors <- c("#4682B4", "tan")


# plot per tree 
hydro_tree_plot <- ggplot(trait_sums_long_hydro, aes(x = Sample_ID, y = CLR_Sum, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = "Relative Over- and Under-Representation of Hydrophobic Traits",
    fill = "Trait"
  ) +
  scale_fill_manual(values=hydro_colors, 
                    name="Trait",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

hydro_tree_plot



## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(trait_sums_long_hydro$Sample_ID)
all_traits <- unique(trait_sums_long_hydro$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
trait_sums_long_hydro_full <- trait_sums_long_hydro %>%
  select(Sample_ID, Trait_clean, CLR_Sum) %>%
  complete(Sample_ID = all_trees, Trait_clean = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 
# Grab site from sample_data
site <- select(sample_data, Sample_ID, Site) 

trait_sums_long_hydro_full <- merge(trait_sums_long_hydro_full, site, by  = "Sample_ID")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

trait_sums_long_hydro_full$Site <- factor(trait_sums_long_hydro_full$Site, levels = site_order)

# and hosts
hosts <- select(sample_data, Sample_ID, Host_ID)

# merge hosts to data 
trait_sums_long_hydro_full <- merge(trait_sums_long_hydro_full, hosts, by = "Sample_ID")


# Clean sample_ID for less visual clutter 
trait_sums_long_hydro_full <- trait_sums_long_hydro_full %>%
  mutate(Tree_ID = Sample_ID %>%
           str_remove("^A-") %>%
           str_remove("^N-") %>%
           str_remove("^S-") %>%
           str_remove("^W-"))

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
hydro_traits_pos <- trait_sums_long_hydro_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(CLR_Sum)) %>%
  ungroup()

hydro_traits_neg <- trait_sums_long_hydro_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup()

# Merge for plotting 

# Combine
hydro_traits_diverging <- bind_rows(hydro_traits_pos, hydro_traits_neg)


# Diverging bar plot for trait representation in individual trees 

hydro_diverging <- ggplot(hydro_traits_diverging, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x") +
  scale_fill_manual(values=hydro_colors, 
                    name="Hyphal Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme_minimal() +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Hydrophobic Traits",
    fill = "Trait_clean") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, color = "black")

hydro_diverging




# Subset to one location and facet by host species 

hydro_north <- dplyr::filter(hydro_traits_diverging, Site == "Northern")


hydro_north_plot <- ggplot(hydro_north, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(values=hydro_colors, 
                    name="Hyphal Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Hydrophobic Traits",
    fill = "Trait_clean") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

hydro_north_plot


hydro_wfdp <- dplyr::filter(hydro_traits_diverging, Site == "WFDP")


hydro_wfdp_plot <- ggplot(hydro_wfdp, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(values=hydro_colors, 
                    name="Hyphal Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Hydrophobic Traits",
    fill = "Trait_clean") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

hydro_wfdp_plot


hydro_andrews <- dplyr::filter(hydro_traits_diverging, Site == "Andrews")


hydro_andrews_plot <- ggplot(hydro_andrews, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(values=hydro_colors, 
                    name="Hyphal Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Hydrophobic Traits",
    fill = "Trait_clean") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"), 
        legend.position = "none")

hydro_andrews_plot


hydro_south <- dplyr::filter(hydro_traits_diverging, Site == "Southern")


hydro_south_plot <- ggplot(hydro_south, aes(x = Tree_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  scale_fill_manual(values=hydro_colors, 
                    name="Hyphal Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- and Under-Representation of Hydrophobic Traits",
    fill = "Trait_clean") +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 11, colour="black"),
    legend.text = element_text(size = 10, colour="black"),
    strip.text = element_text(size = 11, colour="black")) +
  theme(legend.title = element_text(colour="black", size=11, face="bold"), 
        legend.position = "none")

hydro_south_plot


# Remember to run once with legend as "bottom" to get legend for plotting 


# INTERPRETATION:

# These plots are showing the clr-transformed abundance weighted trait values, so not raw abundance 
# A positive value tells us that that trait is relatively over-represented in the community, while 
# a negative value tells us that the trait is relatively under-represented. 


########################################################

# Get weighted taxon abundance in tree communities 

## Using the clr transformed values from the trees_EM_clr object, pulled from the transformed phyloseq

# Keep rownames
tree_ids <- rownames(trees_EM_clr)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_EM_clr <- as.data.frame(trees_EM_clr)
tree_matrix_EM_clr[] <- lapply(tree_matrix_EM_clr, as.numeric)  
tree_matrix_EM_clr <- as.matrix(tree_matrix_EM_clr)
rownames(tree_matrix_EM_clr) <- tree_ids 


# Get dataframe for plotting 
tree_abund_df <- as.data.frame(tree_matrix_EM_clr) %>% rownames_to_column(var = "Sample_code")


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
EM_tax <- tax_table(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

EM_tax <- as.data.frame(EM_tax) %>% rownames_to_column(var = "ASV")

tree_tax_full <- merge(EM_tax, tree_abund_long, by = "ASV")

# Clean up tree names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Tree = Sample_code %>%
           str_replace_all("_FWD", ""))


#### manipulations to plot 

# plot per tree just to get a look at things - not scaled yet 
taxa_tree_plot <- ggplot(tree_tax_full, aes(x = Sample_code, y = CLR_Abund, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "CLR-Weighted Fungal Taxonomic Composition per Tree",
    x = "",
    y = "Relative Genera Abundance",
    fill = "Genus") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 9)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

taxa_tree_plot



## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(tree_tax_full$Sample_code)


# Expand to all Trees and fill missing with 0
# Putting a pin in this, not sure I need to do this 
tree_tax_plotting <- tree_tax_full %>%
  dplyr::select(Sample_code, CLR_Abund, Genus, Family, Species) %>%
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


# Fix family name format 
tree_taxa_diverging <- tree_taxa_diverging %>%
  mutate(Family = Family %>%
           str_remove("f__"))



# Diverging bar plot for taxon representation in individual trees 

tree_fam_plot <- ggplot(tree_taxa_diverging, aes(x = Sample_code, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Tree",
    y = "Relative Over- or Under-Representation of Ectomycorrhizal Families",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 8)) +
  theme_minimal(base_size = 11) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 0, colour="black")) +
  theme(legend.title = element_text(colour="black", size=11, face = "bold"), 
        legend.position = "bottom") +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 11)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 11))


tree_fam_plot



# Subset to one location and facet by host species 

fam_north <- dplyr::filter(tree_taxa_diverging, Site == "Northern")

tree_fam_north_plot <- ggplot(fam_north, aes(x = Sample_code, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Ectomycorrhizal Families",
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

tree_fam_wfdp_plot <- ggplot(fam_wfdp, aes(x = Sample_code, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Ectomycorrhizal Families",
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

tree_fam_andrews_plot <- ggplot(fam_andrews, aes(x = Sample_code, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Ectomycorrhizal Families",
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


tree_fam_south_plot <- ggplot(fam_south, aes(x = Sample_code, y = CLR_Abund_scaled, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Over- or Under-Representation of Ectomycorrhizal Families",
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


## Can return to do genus-level plot if desired 
x


## Done visualizing the functional and taxonomic variation across the sites and hosts 

## -- END -- ## 


