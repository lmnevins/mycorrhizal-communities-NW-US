# -----------------------------------------------------------------------------#
# Differential Abundance Analyses for AM Communities using ALDEx2
# Original Author: L. McKinley Nevins 
# July 14, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     phyloseq v 1.48.0
#                     vegan v 2.6.10
#                     dplyr v 1.1.4
#                     purrr v 1.0.4
#                     ALDEx2 v 1.36.0
#                     ggplot2 v 3.5.1
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
require(tidyverse); packageVersion("tidyverse")
require(phyloseq); packageVersion("phyloseq")
require(vegan); packageVersion("vegan")
require(dplyr); packageVersion("dplyr")
require(purrr); packageVersion("purrr")
require(ALDEx2); packageVersion("ALDEx2")
require(ggplot2); packageVersion("ggplot2")

#################################################################################
#                               Main workflow                                   #
#  Perform differential abundance analyses to determine which taxa are more or  #
#  less abundant at specific sites using ALDEx2.                                #
#                                                                               #
#################################################################################

########
# Load in the environmental and site data 
sample_metadata <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_enviro_all_2025.csv")
# Data for 130 host trees, including all environmental variables along with other grouping variables 

########

# load in AM final phyloseq object 
ps_AM <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_phyloseq_func_subset_final_2025.RDS")

# 363 ASVS in 130 tree samples 

# Pull out ASV table of raw counts from phyloseq object 
asvs <- phyloseq::otu_table(ps_AM) %>% as.data.frame()

asvs <- t(asvs) %>% as.data.frame()


##### 

# Get tax table for later 
tax <- tax_table(ps_AM) %>% as.data.frame() %>% rownames_to_column(var = "ASV")

## Can skip to section 'Environmental Extremes' if only interested in main results in paper 

######Differential Abundance Tests###########################################################################
##Differential abundance using ALDEx2

# Aldex2 performs CLR transformation internally, and this was done each time with a pseudocount of 1, while 
# a pseudocount of 1e-06 was done for the other compositional analyses. This was necessary because aldex.clr 
# applies the clr transformation to each read separately, so it accounts for the pseudocount in a different way. 

# Following some of this vignette: 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html


# ALDEx2 only does comparisons between 2 categories, so to compare between my 4 sites I need 
# to do pairwise comparisons for each pair of sites 


# STEP 1: Ensure ASV matrix is numeric and cleaned
asvs[] <- lapply(asvs, as.numeric)  # convert to numeric
asv_matrix <- as.matrix(asvs)
asv_matrix <- asv_matrix[rowSums(asv_matrix) > 0, ]  # remove all-zero ASVs

# STEP 2: Create named grouping vector
grouping_vector <- sample_metadata$Site
names(grouping_vector) <- sample_metadata$X

# STEP 3: Get all pairwise combinations of sites
site_pairs <- combn(unique(grouping_vector), 2, simplify = FALSE)

# STEP 4: Run ALDEx2 for each pairwise site comparison
pairwise_results <- list()

for (pair in site_pairs) {
  cat("\n----\nRunning ALDEx2 for", pair[1], "vs", pair[2], "\n")
  
  # Find trees from each site
  tree_ids <- names(grouping_vector[grouping_vector %in% pair])
  pair_sites <- grouping_vector[tree_ids]
  
  # Subset ASV matrix to trees from these sites
  sub_asv <- asv_matrix[, tree_ids, drop = FALSE]
  
  # Filter again in case some ASVs are now zero-sum
  sub_asv <- sub_asv[rowSums(sub_asv) > 0, , drop = FALSE]
  
  # Skip if too few samples
  if (length(unique(pair_sites)) < 2 || any(table(pair_sites) < 2)) {
    cat("  ❌ Skipping: not enough samples per site.\n")
    next
  }
  
  # Run ALDEx2
  x.clr <- aldex.clr(sub_asv + 1, conds = pair_sites, mc.samples = 128, denom = "all")
  x.tt <- aldex.ttest(x.clr)
  x.effect <- aldex.effect(x.clr, CI=T)
  
  results <- cbind(x.tt, x.effect)
  results$ASV <- rownames(results)
  results$Comparison <- paste(pair[1], "vs", pair[2])
  
  pairwise_results[[paste(pair[1], pair[2], sep = "_vs_")]] <- results
}

all_results_df <- bind_rows(pairwise_results, .id = "Comparison_ID")


# Merge taxonomic information with the differential abundance results 
aldex_taxa <- all_results_df %>%
  left_join(tax, by = "ASV")



### Exploratory Visualizations ###

# Filter to keep only comparisons with actual data
aldex_taxa <- aldex_taxa %>%
  filter(!is.na(effect)) %>%
  mutate(Significance = ifelse(we.eBH < 0.05 & abs(effect) > 0.6, "Significant", "Not significant"))

# Plot
sig_plot <- ggplot(aldex_taxa, aes(x = diff.btw, y = effect)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey70", "firebrick")) +
  facet_wrap(~ Comparison, scales = "free") +
  theme_minimal() +
  labs(
    title = "Differential Abundance Across Site Comparisons",
    x = "CLR Difference Between Sites",
    y = "Effect Size",
    color = "Significance"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )


sig_plot 

# Assess significant genera

# Summarize number of significant ASVs per genus per comparison
sig_genera <- aldex_taxa %>%
  filter(we.eBH < 0.05, abs(effect) > 0.6) %>%
  group_by(Comparison, Genus) %>%
  summarise(Signif_ASVs = n(), .groups = "drop") %>%
  arrange(desc(Signif_ASVs))

sig_families <- aldex_taxa %>%
  filter(we.eBH < 0.05, abs(effect) > 0.6) %>%
  group_by(Comparison, Family) %>%
  summarise(Signif_ASVs = n(), .groups = "drop") %>%
  arrange(desc(Signif_ASVs))

# ASVs are all in the Glomeraceae


## INTERPRETATION: For these exploratory analyses I was just lumping together all host trees as 
# replicates for each site. We know the communities are different between different 
# host tree taxa, so the comparisons need to take this into account. 

######Analyses considering host tree taxa######################################################

## Analyses considering host tree taxa 

### 1. Comparing host species within each site -> will compare the impact of host 
      # identity while keeping the local environment consistent 

# asvs is the input dataframe where each tree is a column

# sample_metadata identifies each tree to match in the 'X' column 

# Identify hosts by site 
sample_metadata <- sample_metadata %>%
  mutate(Host_Site = paste(Host_ID, Site, sep = "_"))


# Initialize results list
all_results <- list()

# Get unique sites
sites <- unique(sample_metadata$Site)

# Loop over sites
for (s in sites) {
  
  # Subset metadata to current site
  site_meta <- sample_metadata %>% filter(Site == s)
  
  # Get host species at that site
  hosts <- unique(site_meta$Host_ID)
  
  # Get all host species pairs
  host_pairs <- combn(hosts, 2, simplify = FALSE)
  
  # Loop through host pairs
  for (pair in host_pairs) {
    cat("\n----\nRunning ALDEx2 for", pair[1], "vs", pair[2], "\n")
    
    host1 <- pair[1]
    host2 <- pair[2]
    
    # Subset metadata and samples for this pair
    pair_meta <- site_meta %>%
      filter(Host_ID %in% c(host1, host2))
    
    samples <- pair_meta$X
    
    # Subset ASV table
    pair_asv <- asvs[, samples]
    
    # Grouping vector
    group_vector <- pair_meta$Host_ID
    names(group_vector) <- pair_meta$X
    
    # Run ALDEx2
    aldex_result <- aldex.clr(pair_asv +1, conds = group_vector, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex_out <- aldex.ttest(aldex_result, paired.test = FALSE)
    aldex_effect <- aldex.effect(aldex_result, CI = T)
    
    # Combine results
    combined <- cbind(aldex_out, aldex_effect)
    combined$ASV <- rownames(combined)
    combined$Comparison <- paste0(host1, "_vs_", host2)
    combined$Site <- s
    
    # Store
    all_results[[paste(s, host1, host2, sep = "_")]] <- combined
  }
}

# Combine into one dataframe
all_results_host_site <- bind_rows(all_results)


# Merge taxonomic information with the differential abundance results 
all_results_host_site <- all_results_host_site %>%
  left_join(tax, by = "ASV")


# Flag significance 
all_results_host_site <- all_results_host_site %>%
  mutate(Significant = ifelse(we.eBH < 0.05 & abs(effect) > 0.6, "Significant", "Not Significant"),
         Facet_Label = paste0(Site, ": ", Comparison))

# Host x site plot 
host_site_plot <- ggplot(all_results_host_site, aes(x = effect, y = -log10(we.eBH), color = Significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "#D95F02")) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  facet_wrap(~ Facet_Label, scales = "free_y") +
  labs(x = "Effect Size", y = "-log10 Adjusted p-value", color = "Significant",
       title = "Differential Abundance by Host Species Within Sites") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

host_site_plot

## No significant results 


# Summarize number of significant ASVs per genus per comparison
sig_genera2 <- all_results_host_site %>%
  filter(we.eBH < 0.05, abs(effect) > 0.6) %>%
  group_by(Comparison, Genus) %>%
  summarise(Signif_ASVs = n(), .groups = "drop") %>%
  arrange(desc(Signif_ASVs))

sig_families2 <- all_results_host_site %>%
  filter(we.eBH < 0.05, abs(effect) > 0.6) %>%
  group_by(Comparison, Family) %>%
  summarise(Signif_ASVs = n(), .groups = "drop") %>%
  arrange(desc(Signif_ASVs))

# No significant taxa 

######Compare across sites and hosts#######################

### 2. Comparing host species across each site -> will compare how the fungal community of each 
      # host tree changes across the environmental gradients 

# asvs is the input dataframe where each tree is a column

# sample_metadata identifies each tree to match in the 'X' column 

# List of unique host species
host_list <- unique(sample_metadata$Host_ID)
site_pairs_per_host <- list()

# Step 1: Loop through each host species
for (host in host_list) {
  cat("\n----\nRunning ALDEx2 for", pair[1], "vs", pair[2], "\n")
  
  # Subset metadata and ASV table to only this host species
  meta_sub <- sample_metadata %>% filter(Host_ID == host)
  asv_sub <- asvs[, meta_sub$X]
  
  # Get all pairwise site comparisons for this host
  site_pairs <- combn(unique(meta_sub$Site), 2, simplify = FALSE)
  
  for (pair in site_pairs) {
    site1 <- pair[1]
    site2 <- pair[2]
    
    # Subset samples for this site pair
    samples_pair <- meta_sub %>% filter(Site %in% c(site1, site2))
    asv_pair <- asv_sub[, samples_pair$X]
    
    # Grouping vector (site labels)
    group_vector <- samples_pair$Site
    
    # Run ALDEx2
    aldex_result <- aldex.clr(asv_pair +1, group_vector, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex_test <- aldex.ttest(aldex_result)
    aldex_effect <- aldex.effect(aldex_result, CI = T)
    
    # Combine results
    result_df <- cbind(ASV = rownames(aldex_test), aldex_test, aldex_effect) %>%
      mutate(Site_Comparison = paste(site1, "vs", site2, sep = "_"),
             Host_ID = host)
    
    # Save
    site_pairs_per_host[[paste(host, site1, site2, sep = "_")]] <- result_df
  }
}

# Combine all results
all_results_hosts_across <- bind_rows(site_pairs_per_host)

# Merge taxonomic information with the differential abundance results 
all_results_hosts_across <- all_results_hosts_across %>%
  left_join(tax, by = "ASV")


all_results_hosts_across <- all_results_hosts_across %>%
  mutate(Significant = ifelse(we.eBH < 0.05 & abs(effect) > 0.6, "Significant", "Not Significant"))


# Plot host species across the sites 
hosts_across_plot <- ggplot(all_results_hosts_across, aes(x = effect, y = -log10(we.eBH), color = Significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("grey70", "#D7263D")) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.4) +
  facet_grid(Host_ID ~ Site_Comparison, scales = "free_y", switch = "y") +
  labs(
    x = "Effect Size (Centered Log Ratio Difference)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    strip.placement = "outside",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(legend.position = "bottom")

hosts_across_plot

## Dotted lines added to show significance thresholds for both p-value and 
# the effect size. Positive effect size means the taxon is more abundant in the first 
# site, negative means it's more abundant in the second site. 


# Summarize number of significant ASVs per genus per comparison
sig_genera3 <- all_results_hosts_across %>%
  filter(we.eBH < 0.05, abs(effect) > 0.6) %>%
  group_by(Site_Comparison, Host_ID, Genus) %>%
  summarise(Signif_ASVs = n(), .groups = "drop") %>%
  arrange(desc(Signif_ASVs))

sig_families3 <- all_results_hosts_across %>%
  filter(we.eBH < 0.05, abs(effect) > 0.6) %>%
  group_by(Site_Comparison, Host_ID, Family) %>%
  summarise(Signif_ASVs = n(), .groups = "drop") %>%
  arrange(desc(Signif_ASVs))

# All in the Glomeraceae

 
### Conclusions - This second comparison reflects the changes to the fungal communities of each host tree 
# species across the sites


####Environmental Extremes#############################################################################

# Because the environmental gradients aren't necessarily directional across the sites, I also want to 
# find the trees of each species that are at the ends of key environmental gradients, and then assess 
# how the fungal communities for each of those species vary between those ends. 


# Explore the sample_metadata to look at variation in some environmental variables 


## From the forward selection for the dbRDA, the variables selected as significant for the AM taxa were:
# ph + elev + Sand + apr1_SWE


# Soil pH
hist(sample_metadata$ph)

# Elevation 
hist(sample_metadata$elev)

# Soil Sand
hist(sample_metadata$Sand)

# April 1 SWE
hist(sample_metadata$apr1_SWE)

# The range of values will look different for each tree species, but I can target the 
# trees that are at the end of the range for each species independently 


# Test using elevation:
env_var <- "elev"  

# avg env value for each species × site
host_site_env <- sample_metadata %>%
  group_by(Host_ID, Site) %>%
  summarise(mean_env = mean(.data[[env_var]], na.rm = TRUE), .groups = "drop")

# pick site with min and max env values for each species
host_site_extremes <- host_site_env %>%
  group_by(Host_ID) %>%
  filter(mean_env == min(mean_env) | mean_env == max(mean_env)) %>%
  mutate(Env_Group = ifelse(mean_env == min(mean_env), "Low", "High")) %>%
  ungroup()

# merge with tree metadata to get tree-level grouping
tree_group_assignments <- sample_metadata %>%
  inner_join(host_site_extremes, by = c("Host_ID", "Site")) %>%
  dplyr::select(Sample_ID, Host_ID, Site, Env_Group)

# This selects the trees for each species that occur in the high and low sites along their range. 
# Then these can be compared pairwise to see if their communities are different. 


##### 

sample_metadata$sample_code <- sample_metadata$X

# Loop to find the extremes for each species across enviro variables of interest 

# Define environmental variables to analyze
# Looking at main variables related to drought and soil properties that I have a priori
# expectations for being important 

env_vars <- c("elev", "mean_precip_mm", "mean_summer_precip_mm", "MAT", "pct_C", "pct_N", "ph", "org_matter", "Sand", "Silt",
              "Clay", "EC", "tot_bases", "avg_July_SPEI", "count_mod_dry", "count_sev_dry", "apr1_SWE")

# Make an empty list to store results
extreme_tree_groups_list <- list()

# Loop over each environmental variable
for (var in env_vars) {
  
  # Step 1: Calculate mean env per species x site
  avg_env <- sample_metadata %>%
    group_by(Host_ID, Site) %>%
    summarise(mean_env = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  # Step 2: Identify the highest and lowest site for each species
  extremes <- avg_env %>%
    group_by(Host_ID) %>%
    filter(mean_env == max(mean_env) | mean_env == min(mean_env)) %>%
    mutate(Env_Group = ifelse(mean_env == min(mean_env), "Low", "High")) %>%
    ungroup()
  
  # Step 3: Tag the individual trees that match those sites and species
  tree_tags <- sample_metadata %>%
    inner_join(extremes, by = c("Host_ID", "Site")) %>%
    dplyr::select(sample_code, Host_ID, Site, Env_Group) %>%
    mutate(Environmental_Var = var)
  
  # Add this to the list
  extreme_tree_groups_list[[var]] <- tree_tags
}

# Combine everything into a single data frame
extreme_tree_groups <- bind_rows(extreme_tree_groups_list)

###########
# ALDEx2

# asvs is the dataframe of raw sequence reads, with ASVs as rows and tree samples as columns 
# tax is the taxonomic information for each ASV
# column 'sample_code' in the sample_metadata df matches the columns in asvs 

# Initialize results list
aldex_results_list <- list()

# Get unique combinations of Host_ID and Environmental_Var
comparisons <- extreme_tree_groups %>%
  distinct(Host_ID, Environmental_Var)

# Loop through each Host_ID × Environmental_Var pair
for (i in 1:nrow(comparisons)) {
  
  host <- comparisons$Host_ID[i]
  env_var <- comparisons$Environmental_Var[i]
  
  # Subset to relevant samples
  subset_samples <- extreme_tree_groups %>%
    filter(Host_ID == host, Environmental_Var == env_var)
  
  # Make sure both groups are present
  group_counts <- table(subset_samples$Env_Group)
  if (!all(c("High", "Low") %in% names(group_counts)) || any(group_counts < 2)) {
    message(paste("Skipping", host, env_var, "- insufficient group sizes"))
    next
  }
  
  # Extract count data for those samples
  sample_ids <- subset_samples$sample_code
  subset_counts <- asvs[, colnames(asvs) %in% sample_ids]
  
  # Make sure ASVs are rows and Samples are columns
  if (!is.matrix(subset_counts)) {
    subset_counts <- as.matrix(subset_counts)
  }
  
  # Subset and reorder grouping vector to match column order
  group_vector <- subset_samples %>%
    filter(sample_code %in% colnames(subset_counts)) %>%
    arrange(match(sample_code, colnames(subset_counts))) %>%
    pull(Env_Group)
  
  if (length(group_vector) != ncol(subset_counts)) {
    warning(paste("Skipping", host, env_var, "- mismatch in sample alignment"))
    next
  }
  
  # Run ALDEx2
  tryCatch({
    aldex.clr.out <- aldex.clr(subset_counts +1, group_vector, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex.test.out <- aldex.ttest(aldex.clr.out, paired.test = FALSE)
    aldex.effect.out <- aldex.effect(aldex.clr.out, CI=T)
    
    # Combine results
    combined <- cbind(aldex.test.out, aldex.effect.out)
    combined$ASV <- rownames(combined)
    
    # Annotate
    combined <- combined %>%
      mutate(Host_ID = host,
             Environmental_Var = env_var,
             Comparison = paste0("High_vs_Low_", env_var))
    
    # Store
    aldex_results_list[[paste(host, env_var, sep = "_")]] <- combined
    
  }, error = function(e) {
    message(paste("Error for", host, env_var, ":", e$message))
  })
}

# Combine all results 
aldex_all_extremes <- bind_rows(aldex_results_list)


# Merge taxonomic information with the differential abundance results 
aldex_all_extremes <- aldex_all_extremes %>%
  left_join(tax, by = "ASV")


aldex_all_extremes <- aldex_all_extremes %>%
  mutate(Significant = ifelse(we.eBH < 0.05 & abs(effect) > 0.6, "Significant", "Not Significant"))


## Visualizations 

# Plot host species across the sites 
extremes_plot <- ggplot(aldex_all_extremes, aes(x = effect, y = -log10(we.eBH), color = Significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("grey70", "#D7263D")) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.4) +
  facet_grid(Host_ID ~ Environmental_Var, scales = "free_y", switch = "y") +
  labs(
    x = "Effect Size (Centered Log Ratio Difference)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    strip.placement = "outside",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(legend.position = "bottom")

extremes_plot

# Subset the results to just the environmental variables that had a significant result

aldex_all_extremes$Environmental_Var <- as.factor(aldex_all_extremes$Environmental_Var)

signif_comparisons <- aldex_all_extremes %>%
  filter(we.eBH < 0.05 & abs(effect) > 0.6) %>%
  distinct(Environmental_Var, Host_ID)

# Keep all rows from comparisons that had significant hits
aldex_trimmed_for_plot <- aldex_all_extremes %>%
  semi_join(signif_comparisons, by = c("Host_ID", "Environmental_Var"))


# Plot just significant environmental variables 
extremes_plot2 <- ggplot(aldex_trimmed_for_plot, aes(x = effect, y = -log10(we.eBH), color = Significant)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("grey70", "#D7263D")) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.4) +
  facet_grid(Host_ID ~ Environmental_Var, scales = "free_y", switch = "y") +
  labs(
    x = "Effect Size (Centered Log Ratio Difference)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    strip.placement = "outside",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(legend.position = "bottom")

extremes_plot2


# Subset the results table to get what I want for the summary table 

aldex_trimmed_for_plot$Significant <- as.factor(aldex_trimmed_for_plot$Significant)


extremes_sig_results <- filter(aldex_trimmed_for_plot, Significant == "Significant")


# Select just the results columns of interest 
extremes_sig_results <- dplyr::select(extremes_sig_results, we.eBH, effect, ASV, Host_ID, Environmental_Var,
                                      Comparison, Family, Genus)

# Save file 
write.csv(extremes_sig_results, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_extremes_diff_abund.csv")

############################
## Visualize in a more exciting way than just a table 

# Rename a few of the environmental variables so they fit better 
extremes_sig_results$Environmental_Var <- as.factor(extremes_sig_results$Environmental_Var)

# Recode levels
extremes_sig_results <- extremes_sig_results %>%
  mutate(Environmental_Var = dplyr::recode(Environmental_Var,
                                    mean_precip_mm = "MAP",
                                    mean_summer_precip_mm = "MSP",
                                    avg_July_SPEI = "SPEI",
                                    count_mod_dry = "Mod_Dry",
                                    count_sev_dry = "Sev_Dry"))

# Put in alphabetical order for clean plotting 
extremes_sig_results$Environmental_Var <- factor(extremes_sig_results$Environmental_Var, 
                                                 levels=c("SPEI", "Silt", "Sev_Dry", "Sand", 
                                                          "MSP", "Mod_Dry", "MAP", "EC"))

levels(extremes_sig_results$Environmental_Var)


# Try out an ASV label 
extremes_sig_results <- extremes_sig_results %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus),
    ASV_label = paste0(ASV, " (", Family, ", ", Genus, ")")
  )

# This is quite long, so while it may be useful at some point, it's clunky here


# Set ASVs to group by a bit of phylogenetic structure 
extremes_sig_results <- extremes_sig_results %>%
  mutate(ASV = factor(ASV, 
                      levels = extremes_sig_results %>%
                        arrange(Genus, Family, ASV) %>% # Will go by genus first, then family
                        pull(ASV) %>%
                        unique()))


# Find the effect size range to set for plotting 
max_abs_effect <- max(abs(extremes_sig_results$effect), na.rm = TRUE)



###############

# Heat map of effect sizes of changes in abundance 

heat_AM <- ggplot(extremes_sig_results, aes(y = Environmental_Var, 
                   x = ASV, 
                   fill = effect)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, 
                       name = "Effect size") +
  facet_wrap(~ Host_ID, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  
    axis.text.y = element_text(angle = 90, hjust = 1, size = 8),    
    strip.text = element_text(size = 12, face = "bold")) + 
  labs(y = "Environmental gradient", x = "ASV")

heat_AM


# Decent but I don't love it 


## Bubble plot 
bubble_AM <- ggplot(extremes_sig_results, aes(y = Environmental_Var, 
                                              x = ASV, size = abs(effect), color = effect)) +
  geom_point(aes(size = abs(effect), color = effect)) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0,
    limits = c(-3.5, 3.5), name = "Effect size") +
  scale_size_continuous(name = "|Effect size|", limits = c(0, 3.5)) +
  facet_wrap(~ Host_ID, scales = "free_y", nrow = 2) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank()) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  labs(y = "Environmental Factor", x = "")

bubble_AM



## Fine-tuning 


# Subset results to only those environmental variables that were forward selected 

extremes_sig_subset <- filter(extremes_sig_results, Environmental_Var %in% c("pH", "Elev", "Sand"))


# Put in alphabetical order for clean plotting 
extremes_sig_subset$Environmental_Var <- factor(extremes_sig_subset$Environmental_Var, 
                                                 levels=c("Sand", "pH", "Elev"))

levels(extremes_sig_subset$Environmental_Var)


bubble_AM_subset <- ggplot(extremes_sig_subset, aes(y = Environmental_Var, 
                                              x = ASV, size = abs(effect), color = effect)) +
  geom_point(aes(size = abs(effect), color = effect)) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0,
                        limits = c(-3.5, 3.5), name = "Effect size") +
  scale_size_continuous(name = "|Effect size|", limits = c(0, 3.5)) +
  facet_wrap(~ Host_ID, scales = "free_y", nrow = 2) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank()) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  labs(y = "Environmental Factor", x = "")

bubble_AM_subset

# this doesn't feel valuable to do because it removes almost all of the results 


## -- END -- ## 

