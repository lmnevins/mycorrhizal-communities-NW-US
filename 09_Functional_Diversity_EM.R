# -----------------------------------------------------------------------------#
# Chapter 2: Merging EM fungal taxonomic data to functions and explore 
# Original Author: L. McKinley Nevins 
# January 29, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.2.1
#                     fundiversity v 1.1.1
#                     vegan v 2.6.6.1
#                     cluster v 2.1.8
#                     FD v 1.0.12.3
#                     ade4 v 1.7.22
#                     phyloseq v 1.48.0
#                     ape v 5.8.1
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
library(ade4); packageVersion("ade4")
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
#################################################################################
#                               Main workflow                                   #
#  Load in EM fungal taxa table, and the table of functions for genera that     #
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
funcs <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/matched_OTU_funcs_EM.csv")


# this has a categorical 'hydro' column for hydrophobicity, and a 'hydro_binary' column where it's been coded 
# as hydrophilic = 1 and hydrophobic = 0.

# The Exploration Types have all been one-hot encoded as separate binary columns - there were 14 categories


# Pull in phyloseq object with the clr transformed values for all downstream steps 
ps_EM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final.RDS")

# 424 trees and 1,941 ASVs


# reformat funcs a tiny bit just to get OTU's as the species in rownames, and my traits only 
traits_EM <- funcs %>% dplyr::select(OTU2, hydrophilic, hydrophobic, ET_contact, ET_contact_short, ET_short, ET_contact_medium,
                               ET_contact_medium_fringe, ET_contact_medium_smooth, ET_medium_smooth, 
                               ET_medium_fringe, ET_medium_mat, ET_medium_long, ET_medium_long_smooth,
                               ET_medium_long_fringe, ET_contact_long_smooth, ET_long) %>% column_to_rownames(var = "OTU2") 

#OTUs are now the row names and there are columns for the binary coding of each trait level 


################## species by site matrix ########

## Keep the host tree code, because the host tree has it's own community and 
# the functioning likely varies between individual trees. 

# get clr abundance with each host tree as its own 'site' 

#######

trees_EM_clr <- otu_table(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

#this already has each individual tree as it's own row 

# calling this "trees" because each tree is its own site 


####################
## 2. DATA ANALYSIS
# Relative Abundance of Functions 
####################

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
otu_ids <- rownames(traits_EM)

traits_matrix_EM <- as.data.frame(traits_EM)
traits_matrix_EM[] <- lapply(traits_matrix_EM, as.numeric)
traits_matrix_EM <- as.matrix(traits_matrix_EM)
rownames(traits_matrix_EM) <- otu_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_EM)) # None        
sum(is.na(traits_matrix_EM))  # None   

# Make sure OTUs match between traits and trees
# Find shared OTUs between the two datasets 
shared_otus <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM))
# all are shared

# set the OTU's to be in the same order 
otu_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM <- traits_matrix_EM[otu_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM)) #TRUE

# All match, good to go here 

############# 

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

# Can now explore the clr-weighted trait profiles of each host tree! 

# I'm interested in exploration type and hydrophobicity separately, so I want to pull these out of 
# the bigger dataset 

# Exploration types 
exploration_traits <- trait_sums_long %>%
  filter(str_starts(Trait, "ET_"))

hydro_traits <- trait_sums_long %>%
  filter(str_starts(Trait, "hydro"))

# Clean up the trait names 
exploration_traits <- exploration_traits %>%
  mutate(Trait_clean = Trait %>%
           str_remove("^ET_") %>%                # remove prefix
           str_replace_all("_", "-") %>%         # underscores → dashes
           str_to_title()                        # capitalize each word
  )


# Hyphal hydrophobicity
hydro_traits <- hydro_traits %>%
  mutate(Trait_clean = Trait %>%
           str_to_title())



# Set exploration type order to reflect range of short to long distance investment 
exploration_order <- c("Contact", "Contact-Short", "Short", "Contact-Medium", "Contact-Medium-Smooth",
                       "Contact-Medium-Fringe", "Contact-Long-Smooth", "Medium-Smooth", "Medium-Fringe", 
                       "Medium-Mat", "Medium-Long", "Medium-Long-Smooth", "Medium-Long-Fringe", "Long")


# apply order to the trait column 
exploration_traits <- exploration_traits %>%
  mutate(Trait_clean = factor(Trait_clean, levels = exploration_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(14, option = "D", direction = -1)
print(viridis_colors)


# plot per tree 
ET_tree_plot <- ggplot(exploration_traits, aes(x = Sample_ID, y = CLR_Sum, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative Exploration Type Abundance",
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
all_trees <- unique(exploration_traits$Sample_ID)
all_traits <- unique(exploration_traits$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
exploration_traits_full <- exploration_traits %>%
  dplyr::select(Sample_ID, Trait_clean, CLR_Sum) %>%
  complete(Sample_ID = all_trees, Trait_clean = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 
# Grab site from sample_data
site <- dplyr::select(sample_data, Sample_ID, Site) 

exploration_traits_full <- merge(exploration_traits_full, site, by  = "Sample_ID")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

exploration_traits_full$Site <- factor(exploration_traits_full$Site, levels = site_order)

# and hosts
hosts <- dplyr::select(sample_data, Sample_ID, Host_ID)

# merge hosts to data 
exploration_traits_full <- merge(exploration_traits_full, hosts, by = "Sample_ID")



## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos <- exploration_traits_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(CLR_Sum)) %>%
  ungroup()

traits_neg <- exploration_traits_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup() %>%
  mutate(CLR_Sum_scaled = -abs(CLR_Sum_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
traits_diverging <- bind_rows(traits_pos, traits_neg)


# Diverging bar plot for trait representation in individual trees 

ET_diverging <- ggplot(traits_diverging, aes(x = Sample_ID, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_x") +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal() +
  labs(
    x = "Tree",
    y = "Exploration Type Abundance",
    fill = "Exploration Type") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, color = "black")

ET_diverging


#### Traits diverging dataframe has each individual tree, the site and host, and the CLR raw
# and CLR scaled weighted abundance for each level of the exploration type trait 

# Save traits_diverging dataframe to use in future analyses 

write.csv(traits_diverging, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_clr_abund.csv")



#### 
## Calculate averages for each host in each site 

trait_avg <- exploration_traits_full %>%
  group_by(Host_ID, Site, Trait_clean) %>%
  summarise(mean_clr = mean(CLR_Sum, na.rm = TRUE), .groups = "drop")

#Split and scale positive and negative values by host
pos_avg <- trait_avg %>%
  filter(mean_clr > 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(mean_clr)) %>%
  ungroup()

neg_avg <- trait_avg %>%
  filter(mean_clr < 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(abs(mean_clr))) %>%
  ungroup() %>%
  mutate(mean_clr_scaled = -abs(mean_clr_scaled))  # ensure values are negative

trait_avg_scaled <- bind_rows(pos_avg, neg_avg)

# Save traits_avg_scaled dataframe to use in future analyses 
write.csv(trait_avg_scaled, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_avg_scaled.csv")


# PLOT

host_ET_plot <- ggplot(trait_avg_scaled, aes(x = Host_ID, y = mean_clr_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "G", direction = -1), exploration_order)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme_minimal() +
  labs(
    y = "Relative Exploration Type Abundance",
    fill = "Exploration Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
theme(legend.title = element_text(colour="black", size=12, face="bold"))

host_ET_plot



######## hydrophobicity analyses ######

# set hydro palette
# hydrophilic   hydrophobic
hydro_colors <- c("#4682B4", "tan")


# plot per tree 
hydro_tree_plot <- ggplot(hydro_traits, aes(x = Sample_ID, y = CLR_Sum, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "CLR-Weighted Hydrophilic Trait Composition per Tree",
    x = "",
    y = "Relative Trait Abundance",
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
all_trees <- unique(hydro_traits$Sample_ID)
all_traits <- unique(hydro_traits$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
hydro_traits_full <- hydro_traits %>%
  select(Sample_ID, Trait_clean, CLR_Sum) %>%
  complete(Sample_ID = all_trees, Trait_clean = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 
# Grab site from sample_data
site <- select(sample_data, Sample_ID, Site) 

hydro_traits_full <- merge(hydro_traits_full, site, by  = "Sample_ID")

site_order <- c("Northern", "WFDP", "Andrews", "Southern")

hydro_traits_full$Site <- factor(hydro_traits_full$Site, levels = site_order)

# and hosts
hosts <- select(sample_data, Sample_ID, Host_ID)

# merge hosts to data 
hydro_traits_full <- merge(hydro_traits_full, hosts, by = "Sample_ID")



## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
hydro_traits_pos <- hydro_traits_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(CLR_Sum)) %>%
  ungroup()

hydro_traits_neg <- hydro_traits_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Sample_ID) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup() %>%
  mutate(CLR_Sum_scaled = -abs(CLR_Sum_scaled))  # ensure values are negative

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
    x = "Tree",
    y = "Relative Hydrophobicity Trait Abundance",
    fill = "Trait_clean") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, color = "black")

hydro_diverging


## Calculate averages for each host in each site 

hydro_trait_avg <- hydro_traits_full %>%
  group_by(Host_ID, Site, Trait_clean) %>%
  summarise(mean_clr = mean(CLR_Sum, na.rm = TRUE), .groups = "drop")

#Split and scale positive and negative values by host
hydro_pos_avg <- hydro_trait_avg %>%
  filter(mean_clr > 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(mean_clr)) %>%
  ungroup()

hydro_neg_avg <- hydro_trait_avg %>%
  filter(mean_clr < 0) %>%
  group_by(Site, Host_ID) %>%
  mutate(mean_clr_scaled = mean_clr / sum(abs(mean_clr))) %>%
  ungroup() %>%
  mutate(mean_clr_scaled = -abs(mean_clr_scaled))  # ensure values are negative

hydro_trait_avg_scaled <- bind_rows(hydro_pos_avg, hydro_neg_avg)


# PLOT

host_hydro_plot <- ggplot(hydro_trait_avg_scaled, aes(x = Host_ID, y = mean_clr_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  scale_fill_manual(values=hydro_colors, 
                    name="Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "Tree",
    y = "Relative Hydrophobicity Trait Abundance",
    fill = "Trait_clean") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    legend.text = element_text(size = 11, colour="black"),
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face="bold"))

host_hydro_plot


# INTERPRETATION:

# These plots are showing the clr-transformed abundance weighted trait values, so not raw abundance 
# A positive value tells us that that trait is relatively more abundant than most of the other OTU's
# in the community, while a negative value tells us that the trait is relatively less abundant, 
# because the OTUs that have that trait are relatively less abundant than other OTU's. 


########################################################

# Get weighted taxon abundance in tree communities 

# Get dataframe for plotting 
tree_abund_df <- as.data.frame(tree_matrix_EM) %>% rownames_to_column(var = "Sample_code")


# Join table of tree environmental data 

#subset sample data 
sample_data <- dplyr::select(sample_data, Sample_code, Site, Host_ID)


tree_abund_full <- tree_abund_df %>%
  merge(sample_data, by = "Sample_code")


# Reshape to long format
tree_abund_long <- tree_abund_full %>%
  pivot_longer(-c(Sample_code, Site, Host_ID), names_to = "OTU", values_to = "CLR_Abund")


## Merge with taxa data 
# Pull out tax table from the trimmed table for taxa that have functional assignments 
EM_tax <- tax_table(ps_EM_clr) %>% as("matrix") %>% as.data.frame()

EM_tax <- as.data.frame(EM_tax) %>% rownames_to_column(var = "OTU")

tree_tax_full <- merge(EM_tax, tree_abund_long, by = "OTU")


# Clean up taxa names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Genus = Genus %>%
           str_replace_all("g__", "")) %>%
  mutate(Species = Species %>%
           str_replace_all("s__", ""))

# Clean up tree names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Tree = Sample_code %>%
           str_replace_all("_FWD", ""))


#### manipulations to plot 

# plot per tree just to get a look at things 
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
  ungroup() %>%
  mutate(CLR_Abund_scaled = -abs(CLR_Abund_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
tree_taxa_diverging <- bind_rows(tree_taxa_pos, tree_taxa_neg)


# Diverging bar plot for taxon representation in individual trees 

tree_taxa_diverging <- ggplot(tree_taxa_diverging, aes(x = Sample_code, y = CLR_Abund_scaled, fill = Genus)) +
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
  group_by(Host_ID, Site, Genus) %>%
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

# wrap some long genera in the legend 
tree_taxa_avg_scaled <- tree_taxa_avg_scaled %>%
  dplyr::mutate(Genus = stringr::str_wrap(Genus, width = 15))

host_tree_taxa_plot <- ggplot(tree_taxa_avg_scaled, aes(x = Host_ID, y = mean_clr_scaled, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  theme_minimal() +
  labs(y = "Relative Genus Abundance",
    fill = "Genus") +
  theme_minimal(base_size = 12) +
  guides(fill = guide_legend(ncol = 4, label.theme = element_text(lineheight = 0.8))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 12)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1)

host_tree_taxa_plot



## Get a summary table of just the relatively most abundant genera for each host and site 

clr_summary <- tree_taxa_avg_scaled %>%
  filter(mean_clr > 0) %>%                             
  group_by(Site, Host_ID, Genus) %>%
  summarise(mean_CLR = mean(mean_clr, na.rm = TRUE)) %>% 
  arrange(Site, Host_ID, desc(mean_CLR)) %>%
  group_by(Site, Host_ID) %>%
  ungroup()


summary_plot_genera  <- ggplot(clr_summary, aes(x = reorder(Genus, mean_CLR), y = mean_CLR, fill = Genus)) +
  geom_col(show.legend = FALSE) +
  facet_grid(Site ~ Host_ID, scales = "free_y") +
  labs(y = "Mean positive CLR", x = "Genus") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
  theme_minimal()

summary_plot_genera



# average for families instead 

# Clean up names a bit 
tree_tax_plotting <- tree_tax_plotting %>%
  mutate(Family = Family %>%
           str_replace_all("f__", ""))


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
  labs(y = "Relative Family Abundance",
    fill = "Family") +
  guides(fill = guide_legend(ncol = 2)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, colour="black")) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(strip.text = element_text(size = 11, colour="black")) +
  theme(axis.title = element_text(colour="black", size = 12)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1)

host_tree_taxa_fam_plot



########################################################################################

#####################
## 3. DATA ANALYSIS
# PERMANOVA for traits x host 
# Beta dispersion of functions
##########@@########

# Load in the position of each host tree on the two PCA axes of the environmental variables 
tree_scores <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_tree_scores_PCs.csv")

tree_scores$Sample_ID <- gsub("_FWD", "", tree_scores$Tree)

# Tree scores has the position of each tree on the two axes, but because all trees in a site are represented 
# by the same environmental information, these scores can be summarized to the tree x site level. 

# Add columns for Host_ID and Site to this dataframe 

dat <- dplyr::select(sample_data, Host_ID, Site, Tree = Sample_code)

tree_scores <- merge(tree_scores, dat, by = "Tree")

tree_scores <- tree_scores %>%
  group_by(Site, Host_ID) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

# Create a site-Host ID column for merging 
tree_scores <- tree_scores %>% mutate(id_code = paste(Site, Host_ID, sep = "_"))

# Trim Host_ID column so they are not duplicated when merging 
tree_scores <- dplyr::select(tree_scores, PC1, PC2, id_code)


# Grab the scaled average CLR-weighted abundance values of each exploration type for each host species across sites 
mod_traits <- trait_avg_scaled

# Create a site-Host ID column for merging 
mod_traits <- mod_traits %>% mutate(id_code = paste(Site, Host_ID, sep = "_"))

# Trim site column so they are not duplicated when merging 
mod_traits <- dplyr::select(mod_traits, Trait_clean, mean_clr, mean_clr_scaled, id_code, Host_ID)

# Merge dataframes by id_code column 
mod_data <- merge(tree_scores, mod_traits, by = "id_code")


## Dataset now has the position of each tree along the two PCA axes of variation for the environmental variables, 
# and the mean CLR scaled abundance for each of the exploration types for each host species in each site 


# Also make a dataset for each individual tree instead 
indiv_traits <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_clr_abund.csv")

# reformat to get abundance for each trait and the sample_ID as the rownames 

  indiv_matrix <- indiv_traits %>%
    mutate(CLR_Sum_scaled = replace_na(CLR_Sum_scaled, 0)) %>%
    pivot_wider(
      id_cols = Sample_ID,
      names_from = Trait_clean,
      values_from = CLR_Sum_scaled,
      values_fill = 0  
    ) %>%
    column_to_rownames("Sample_ID")

## now this dataset is a trait matrix of the CLR weighted abundance values for each exploration type present 
  # in each tree community. Now I can analyze the variation at the individual tree level and how the functional 
  # variation clusters 


# save indiv_matrix
  
write.csv(indiv_matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_tree_ET_clr.csv")
  
  
## Get data summarized at the site_host level 

func_matrix <- mod_data %>%
  tidyr::pivot_wider(
    id_cols = id_code,
    names_from = Trait_clean,
    values_from = mean_clr_scaled,
    values_fill = 0  # fill missing types with 0 (absence)
  ) %>%
  tibble::column_to_rownames("id_code")  # set rownames as id_code


# Create response matrix (just the exploration types)
trait_mat <- func_matrix %>%
  dplyr::select(`Medium-Long`, `Medium-Long-Fringe`, `Contact-Short`, Short, `Contact-Medium`, Contact, 
`Contact-Medium-Smooth`, `Medium-Smooth`, `Medium-Fringe`, `Medium-Long-Smooth`, Long, `Contact-Long-Smooth`,
`Medium-Mat`, `Contact-Medium-Fringe`) %>%
  as.matrix()


# Create metadata
meta <- mod_data %>% 
  dplyr::select(Site, Host_ID, id_code) 


meta <- meta %>% distinct() %>% column_to_rownames("id_code") 

meta_data <- rownames_to_column(meta, var = "id_code")


# Get metadata for the individual tree level 

indiv_meta <- indiv_traits %>% 
  dplyr::select(Site, Host_ID, Sample_ID) 


indiv_meta <- indiv_meta %>% distinct() %>% column_to_rownames("Sample_ID") 

indiv_meta_data <- rownames_to_column(indiv_meta, var = "Sample_ID")


# calculate Aitchison distance using dist() from base R 

aitchison_EM_traits <- dist(indiv_matrix, method = "euclidean")

# Beta diversity of each tree in functional space 

# save Aitchison Distance matrix 
save(aitchison_EM_traits, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_traits.Rdata")


write.csv(mat, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_traits_check.csv")

## visualize ##

# PCA is appropriate for euclidean distances 
# Use the indiv_matrix with clr weighted trait values 

#### Perform PCA on the exploration type traits for each individual tree 
traits_for_pca_EM <- merge(indiv_meta, indiv_matrix, by = 'row.names')


#PCA of differences in composition for Sites and Hosts
pca_traits_EM = prcomp(traits_for_pca_EM[4:17], center = T, scale = F)

sd.pca_traits_EM = pca_traits_EM$sdev
loadings.pca_traits_EM = pca_traits_EM$rotation
names.pca_traits_EM = colnames(traits_for_pca_EM[4:17])
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
                   # ABAM        ABGR        ABPR          ALRU         PSME         TABR         THPL       TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66")

sites <- c(15,16,17,18)


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
                      breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "THPL", "TSHE")) +
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

PCA_traits_site_EM

# Plot the Results by host alone
PCA_traits_host_EM <- ggplot(scores.pca_traits_EM, aes(x = PC1, y = PC2, color = Host_ID)) +
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

PCA_traits_host_EM


#####assess multivariate homogeneity of sites###########
#betadisper function in vegan 

#first object needs to be a dist object of Aitchison distances 
#second object is the groups of interest, as a vector
betadisper_traits_site <- vegan::betadisper(aitchison_EM_traits, indiv_meta_data$Site, type = "median", sqrt.dist = FALSE)

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
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
# Groups      3  1.214 0.40474 4.6758    999  0.003 **
#   Residuals 408 35.317 0.08656                        
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
# diff         lwr       upr     p adj
# Northern-Andrews   0.04421403 -0.05722276 0.1456508 0.6746641
# Southern-Andrews   0.13041232  0.02497869 0.2358460 0.0082901
# WFDP-Andrews       0.12038894  0.01433873 0.2264392 0.0187968
# Southern-Northern  0.08619829 -0.02010312 0.1924997 0.1574790
# WFDP-Northern      0.07617491 -0.03073808 0.1830879 0.2570373
# WFDP-Southern     -0.01002338 -0.12073569 0.1006889 0.9955093

# Significant differences between 
# Southern - Andrews
# WFDP - Andrews


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

centroid_plot_EM


########

#check betadispersion between host tree taxa 
betadisper_traits_host <- betadisper(aitchison_EM_traits, indiv_meta_data$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper_traits_host

permutest(betadisper_traits_host)
#groups are different

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df Sum Sq  Mean Sq     F N.Perm Pr(>F)   
# Groups      7  1.903 0.271788 3.116    999  0.005 **
#   Residuals 404 35.238 0.087223                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 


#if the dispersion is different between groups, then examine
plot(betadisper_traits_host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.8, col = c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66"))

boxplot(betadisper_traits_host)
mod.HSD_traits_host <- TukeyHSD(betadisper_traits_host)
mod.HSD_traits_host
plot(mod.HSD_traits_host)

# $group
#                 diff         lwr          upr     p adj
# ABGR-ABAM -0.023526257 -0.19218871  0.145136202 0.9998847
# ABPR-ABAM -0.076304955 -0.25673681  0.104126895 0.9028293
# ALRU-ABAM -0.001796424 -0.18610050  0.182507655 1.0000000
# PSME-ABAM  0.035887776 -0.13050186  0.202277415 0.9979695
# TABR-ABAM -0.140886502 -0.31481117  0.033038166 0.2122576
# THPL-ABAM -0.171371434 -0.34434335  0.001600480 0.0541950
# TSHE-ABAM -0.073323263 -0.23971290  0.093066376 0.8819047
# ABPR-ABGR -0.052778699 -0.23595565  0.130398248 0.9878822
# ALRU-ABGR  0.021729833 -0.16526251  0.208722171 0.9999666
# PSME-ABGR  0.059414032 -0.10994846  0.228776523 0.9628447
# TABR-ABGR -0.117360245 -0.29413109  0.059410599 0.4676789
# THPL-ABGR -0.147845178 -0.32367869  0.027988336 0.1731906
# TSHE-ABGR -0.049797007 -0.21915950  0.119565484 0.9863444
# ALRU-ABPR  0.074508531 -0.12316484  0.272181905 0.9455083
# PSME-ABPR  0.112192731 -0.06889366  0.293279122 0.5601855
# TABR-ABPR -0.064581547 -0.25261495  0.123451853 0.9669090
# THPL-ABPR -0.095066479 -0.28221896  0.092086005 0.7808813
# TSHE-ABPR  0.002981692 -0.17810470  0.184068083 1.0000000
# PSME-ALRU  0.037684200 -0.14726072  0.222629115 0.9985900
# TABR-ALRU -0.139090078 -0.33084226  0.052662106 0.3478171
# THPL-ALRU -0.169575010 -0.36046344  0.021313421 0.1235027
# TSHE-ALRU -0.071526839 -0.25647176  0.113418076 0.9377433
# TABR-PSME -0.176774278 -0.35137788 -0.002170674 0.0448008
# THPL-PSME -0.207259210 -0.38091379 -0.033604635 0.0074730
# TSHE-PSME -0.109211039 -0.27631023  0.057888154 0.4889991
# THPL-TABR -0.030484932 -0.21137216  0.150402297 0.9995922
# TSHE-TABR  0.067563239 -0.10704037  0.242166842 0.9375732
# TSHE-THPL  0.098048171 -0.07560640  0.271702746 0.6739890

# significant differences are 

# TABR-PSME
# THPL-PSME


# Nicer boxplot 
distances_EM2 <- data.frame(
  Host = betadisper_traits_host$group,
  DistanceToCentroid = betadisper_traits_host$distances
)

centroid_plot_EM2 <- ggplot(distances_EM2, aes(x = Host, y = DistanceToCentroid, fill = Host)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal() +
  labs(x = "Host",
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


# Run PERMANOVA models - for the individual tree level effects of both host and site 
# ** This is permutational so the results will be a little different every time. 

permanova.site <- vegan::adonis2(indiv_matrix ~ Site, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.site

# significant, and only 2.1% explained


permanova.host <- adonis2(indiv_matrix ~ Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.host

# significant, and only 2.7% explained

permanova.both <- vegan::adonis2(indiv_matrix ~ Site * Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both

# significant, and only 11.5% explained


permanova.both2 <- vegan::adonis2(indiv_matrix ~ Site + Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both2

# significant, and only 4.9% explained


########################################################################################
