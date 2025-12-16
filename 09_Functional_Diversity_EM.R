# -----------------------------------------------------------------------------#
# Chapter 2: Merging EM fungal taxonomic data to functions and explore 
# Original Author: L. McKinley Nevins 
# January 29, 2025
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
#  analyses of functional variation across sites and host species, including    #
#  7 EMF hosts: ABAM, ABGR, ABPR, ALRU, PSME, TABR, TSHE                        #
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


# Pull in phyloseq object with the clr transformed values for all downstream steps 
ps_EM_clr <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_transformed_final_no_THPL.RDS")

# 371 trees and 1,627 ASVs


# reformat funcs a tiny bit just to get ASVs as the species in rownames, and my traits only 
traits_EM <- funcs %>% dplyr::select(ASV2, hydrophilic, hydrophobic, ET_contact, ET_contact_short, ET_short, ET_contact_medium,
                               ET_contact_medium_fringe, ET_contact_medium_smooth, ET_medium_smooth, 
                               ET_medium_fringe, ET_medium_mat, ET_medium_long, ET_medium_long_smooth,
                               ET_medium_long_fringe, ET_contact_long_smooth, ET_long) %>% column_to_rownames(var = "ASV2") 

#ASVs are now the row names and there are columns for the binary coding of each trait level 


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

write.csv(traits_diverging, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_clr_abund_no_THPL.csv")



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
write.csv(trait_avg_scaled, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_avg_scaled_no_THPL.csv")


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
    x = "",
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
# A positive value tells us that that trait is relatively more abundant than most of the other ASVs
# in the community, while a negative value tells us that the trait is relatively less abundant, 
# because the ASVs that have that trait are relatively less abundant than other ASVs. 


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
indiv_traits <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_ET_clr_abund_no_THPL.csv")

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
  
write.csv(indiv_matrix, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_tree_ET_clr_no_THPL.csv")
  
  
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
save(aitchison_EM_traits, file="~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_traits_no_THPL.Rdata")

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
                   # ABAM        ABGR        ABPR          ALRU         PSME         TABR     TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#8B0AA5FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#fffd66")

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
#                  Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
#      Groups      3  1.2701 0.42336 5.0857    999  0.002 **
#   Residuals 358 29.8019 0.08325                        
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
#                          diff         lwr       upr     p adj
# Northern-Andrews   0.04219161 -0.064202206 0.14858542 0.7357404
# Southern-Andrews   0.14923028  0.039108390 0.25935217 0.0029566
# WFDP-Andrews       0.11696768  0.005392412 0.22854296 0.0357839
# Southern-Northern  0.10703867 -0.003606425 0.21768377 0.0620696
# WFDP-Northern      0.07477607 -0.037315619 0.18686777 0.3137894
# WFDP-Southern     -0.03226260 -0.147898815 0.08337362 0.8890461

# Significant differences between 
# Southern - Andrews = 0.003 *
# WFDP - Andrews = 0.036 *


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
  labs(x = "",
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
#                  Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#      Groups      6  1.2862 0.214374 2.4605    999  0.017 *
#   Residuals 355 30.9303 0.087128                          
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

# $group
#                 diff         lwr          upr     p adj
# ABGR-ABAM -0.011771788 -0.17582955  0.15228598 0.9999922
# ABPR-ABAM -0.073588216 -0.24909405  0.10191762 0.8766310
# ALRU-ABAM  0.005990242 -0.17328211  0.18526259 0.9999999
# PSME-ABAM  0.030642930 -0.13120406  0.19248992 0.9977832
# TABR-ABAM -0.149723983 -0.31890029  0.01945233 0.1217811
# TSHE-ABAM -0.088364622 -0.25021162  0.07348237 0.6701413
# ABPR-ABGR -0.061816427 -0.23999242  0.11635956 0.9471235
# ALRU-ABGR  0.017762031 -0.16412518  0.19964925 0.9999516
# PSME-ABGR  0.042414718 -0.12232397  0.20715340 0.9881654
# TABR-ABGR -0.137952195 -0.30989698  0.03399259 0.2105337
# TSHE-ABGR -0.076592833 -0.24133152  0.08814585 0.8129631
# ALRU-ABPR  0.079578458 -0.11269819  0.27185510 0.8832187
# PSME-ABPR  0.104231146 -0.07191136  0.28037365 0.5796120
# TABR-ABPR -0.076135768 -0.25903562  0.10676409 0.8803319
# TSHE-ABPR -0.014776406 -0.19091891  0.16136610 0.9999803
# PSME-ALRU  0.024652688 -0.15524300  0.20454838 0.9996490
# TABR-ALRU -0.155714226 -0.34223134  0.03080289 0.1713778
# TSHE-ALRU -0.094354864 -0.27425055  0.08554083 0.7108618
# TABR-PSME -0.180366913 -0.35020362 -0.01053020 0.0291749
# TSHE-PSME -0.119007552 -0.28154473  0.04352963 0.3139174
# TSHE-TABR  0.061359361 -0.10847735  0.23119607 0.9361063

# significant differences are 

# TABR-PSME = 0.029


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
                    breaks=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE"),
                    labels=c("ABAM", "ABGR", "ABPR", "ALRU", "PSME", "TABR", "TSHE")) +
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

# significant, and only 2.2% explained


permanova.host <- adonis2(indiv_matrix ~ Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.host

# not significant, and only 2.4% explained

permanova.both <- vegan::adonis2(indiv_matrix ~ Site * Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both

# significant, and only 10.8% explained


permanova.both2 <- vegan::adonis2(indiv_matrix ~ Site + Host_ID, data = indiv_meta, method = "euclidean", permutations = 999)
permanova.both2

# significant, and only 4.7% explained


########################################################################################


## -- END -- ## 
