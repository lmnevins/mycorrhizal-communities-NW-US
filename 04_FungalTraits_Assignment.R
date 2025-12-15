# -----------------------------------------------------------------------------#
# Use FungalTraits to assign functional guilds to ITS ASVs
# Original Author: L. McKinley Nevins 
# August 7, 2025 
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     phyloseq v 1.48.0
#                     vegan v 2.6.10
#                     fungaltraits v 0.0.3
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####

install.packages("devtools")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(fungaltraits); packageVersion("fungaltraits")


#################################################################################
#                               Main workflow                                   #
#  Use the FungalTraits database through R to assign guilds to ITS ASVs. For    #
#  ASVs that do not have assignments, run through the FunGuild pipeline to see  #
#  if we can get better coverage for assignments.                               #
#                                                                               #
#################################################################################


################ -- 
# (1) DATA PREP
################ -- 

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/"
setwd(wd)


# Load in ITS data that has been cleaned to remove non-fungi, but hasn't undergone other steps 
# this is an OTU table 

ITS_otus <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FUNGuild/EM_OTU_table.csv") 
# This has 17,879 ASVs that were identified, and right now rownames are the sequences, and the last column 
# 'taxonomy' contains all of the taxonomic assignments as a single column separated by semi-colons 


# Need to restructure the taxonomy column to get all separate taxonomic ranks 

ITS_otus <- ITS_otus %>%
  separate(
    col = taxonomy,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    fill = "right",   
    remove = FALSE    
  )


# Load FungalTraits data into R for analyses 

fungal_traits()
traits <- fungal_traits()

# FungalTraits has a Genus column, and since guilds are generally conserved at the genus level, I want to 
# make assignments here. So I need to reformat my genus names to get them to match 

ITS_otus$Genus <- sub("^g__", "", ITS_otus$Genus)
  
# Count unique genera in the otu table 
gen_count_otu <- ITS_otus %>%
  dplyr::count(Genus) %>%
  filter(n > 1)
# 401 genera, some appearing many times and others only once 
# 6,827 have no genus assigned 

# Count species assigned 
spec_count_otu <- ITS_otus %>%
  dplyr::count(Species) %>%
  filter(n > 1)
# 377 species, majority don't have species assignments, but for the ones that do I should match traits 
# according to species for max accuracy 

# Reformat species column for matching with the 'species' column in traits2, which actually 
# has the genus_species listed 
ITS_otus$Species <- sub("^s__", "", ITS_otus$Species)


# add lowercase 'species' column with genus and species together 
ITS_otus <- ITS_otus %>%
  mutate(species = paste(Genus, Species, sep = "_"))



# Trim the traits columns down to just those of interest 
traits2 <- dplyr::select(traits, species, speciesMatched, ifungorum_number, confidence_fg, Genus, growth_form_fg, 
                         guild_fg, higher_clade, substrate, taxonomic_level_fg, trait_fg, trophic_mode_fg, 
                         )
  
  
# Explore how many trait entries exist for each genus 

gen_count_traits <- traits2 %>%
  dplyr::count(Genus) %>%
  filter(n > 1)

# 493 unique genera, ones with the most trait entries are Cladonia, Cortinarius, and Russula


# All of the duplicates seem to be for data in other columns that I already trimmed, so I should 
# be able to summarize this down to just unique values in the 'species' column, which will leave 
# duplicates in the Genus column, but they will be easier to sort through 

# Keep just unique 'species'
Traits_unique <- traits2 %>%
  distinct(species, .keep_all = TRUE)


### Data prep done - the Traits_unique dataframe has 3,591 unique Genus_Species IDs in the 'species' column, 
# and the ITS_otus has the right columns for matching 

#################################################################################

###################### -- 
# (2) GUILD MATCHING 
###################### -- 

# Merge ITS_otus and traits2 to get the functional assignments
# Many ASVs were not assigned to the genus level, so a big fraction won't get assignments here, but 
# we have to go to the genus level to be accurate 


# 1. Try to match ASVs using just the Genus
otu_with_traits_genus <- ITS_otus %>%
  left_join(Traits_unique, by = "Genus")

# NOTE: many-to-many join is expected at the genus level


# 2. Keep only rows with guild assignments
otu_with_traits_genus$guild_fg <- as.factor(otu_with_traits_genus$guild_fg)

traits_genus_trim <- otu_with_traits_genus %>%
  filter(!is.na(guild_fg))
# ~10,896 rows with genus-level guild assignments


# 3. Restrict to ectomycorrhizal-related guilds
# There are 6 classifications that include 'Ectomycorrhizal' in the name 
traits_genus_ecto <- traits_genus_trim %>%
  filter(guild_fg %in% c(
    "Ectomycorrhizal",
    "Ectomycorrhizal-Endophyte-Ericoid Mycorrhizal-Litter Saprotroph-Orchid Mycorrhizal",
    "Ectomycorrhizal-Fungal Parasite",
    "Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph",
    "Ectomycorrhizal-Undefined Saprotroph",
    "Ectomycorrhizal-Wood Saprotroph"
  ))
# ~5,432 rows (still duplicated by OTU)


# 4. An OTU is considered valid if it has AT LEAST ONE
# non-Possible ectomycorrhizal assignment

otu_validity <- traits_genus_ecto %>%
  group_by(OTU) %>%
  summarise(
    has_valid_assignment = any(confidence_fg != "Possible"),
    .groups = "drop"
  )

valid_otus <- otu_validity %>%
  filter(has_valid_assignment) %>%
  pull(OTU)


# 5. Keep only valid OTUs, then drop low-confidence rows
traits_genus_ecto_valid <- traits_genus_ecto %>%
  filter(OTU %in% valid_otus) %>%
  filter(confidence_fg != "Possible")


# 6. collapse to one row per OTU
final_otu_traits <- traits_genus_ecto_valid %>%
  distinct(OTU, .keep_all = TRUE)

# ~2,140 retained


# Checks
final_otu_traits$Genus <- as.factor(final_otu_traits$Genus)
final_otu_traits$confidence_fg <- as.factor(final_otu_traits$confidence_fg)

summary(final_otu_traits$confidence_fg)


# No low-confidence assignments retained
final_otu_traits %>%
  filter(confidence_fg == "Possible") %>%
  nrow()

# No NA guilds
final_otu_traits %>%
  filter(is.na(guild_fg)) %>%
  nrow()


# Save dataset for trait assignments from the literature 
write.csv(final_otu_traits, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/EM_FungalTraits_final.csv")


## Conclusion: 2,053 ASVs were given guild assignments to be ectomycorrhizal. Can subset the ones that have been ID'ed 
# to the genus level but did not have functional guild assignments made, and run these through FunGuild to see if we 
# can get more coverage. 


#################################################################################

########################### -- 
# (3) SUBSET FOR FUNGUILD 
########################### -- 

# Subset ones with genus info but no assignment, or which had 'possible' as the 
# confidence level so they can be checked with FUNGuild instead.

# In the otu_with_traits_genus dataframe, only ASVs with an assignment are given a confidence_fg rank, 
# and ones with no rank just have 'NA'. Therefore, this can be filtered to select ASVs with either 
# 'Possible' or 'NA' ranks to be processed through FUNGuild 

# Use 'otu_with_traits_genus' df because this one shows which ASVs have genus info, and then 
# can be filtered for which don't have guild assignments 

# Check assignment status for all ASVs
otu_assignment_status <- otu_with_traits_genus %>%
  group_by(OTU) %>%
  summarise(
    has_high_confidence = any(confidence_fg %in% c("Probable", "Highly Probable")),
    .groups = "drop"
  ) %>%
  mutate(
    to_reprocess = !has_high_confidence
  )

summary(otu_assignment_status)


# find ASVs that need reprocessing
reprocess_otus <- otu_assignment_status %>%
  filter(to_reprocess) %>%
  pull(OTU)


# create final no guild dataset 
final_no_guild <- otu_with_traits_genus %>%
  filter(OTU %in% reprocess_otus) %>%
  distinct(OTU, .keep_all = TRUE)

# Checks 

# No high-confidence assignments in reprocess set
otu_with_traits_genus %>%
  filter(OTU %in% final_no_guild$OTU,
         confidence_fg %in% c("Probable", "Highly Probable")) %>%
  nrow()

# No overlap with assigned set
intersect(final_no_guild$OTU, final_otu_traits$OTU)

# 13,458 ASVs


## Needs to be in particular format for FUNGuild, can find info about it here: 
# https://github.com/UMNFuN/FUNGuild

# The df already has the taxonomy column FUNGuild needs, and has OTUs in rows and the samples in columns 

# Just need to remove the excess columns and save 

final_no_guild <- dplyr::select(final_no_guild, -c(Kingdom, Phylum, Class, Order, Family, 
                                                               Genus, Species, species.x, species.y, speciesMatched,
                                                               ifungorum_number, confidence_fg, growth_form_fg, guild_fg,
                                                               higher_clade, substrate, taxonomic_level_fg, trait_fg, 
                                                               trophic_mode_fg))
#save dataframe as a csv file 
write.csv(final_no_guild, "./FUNGuild/ITS_OTU_table_2025.csv")

#saves with rownumbers in the first column. After it saves, open it up in Excel just 
#to doublecheck it, and delete that first column. 

## -- END -- ## 
