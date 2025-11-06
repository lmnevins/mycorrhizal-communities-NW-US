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
#    Use the FungalTraits database through R to assign guilds to ITS ASVs and   #
#    then check how this compares to the previous assignments made by FungGuild #
#                                                                               #
#################################################################################


###############
# (1) DATA PREP
###############

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
  count(Genus) %>%
  filter(n > 1)
# 401 genera, some appearing many times and others only once 
# 6,827 have no genus assigned 

# Count species assigned 
spec_count_otu <- ITS_otus %>%
  count(Species) %>%
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
  count(Genus) %>%
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

###################################################################################################

# Merge ITS_otus and traits2 to get the functional assignments
# Many ASVs were not assigned to the genus level, so a big fraction won't get assignments here, but 
# we have to go to the genus level to be accurate 


# For accuracy, need to do a couple rounds here 

# 1. First, match on full genus_species name to be most accurate, using 'species' column in both 
otu_with_traits_species <- ITS_otus %>%
  left_join(Traits_unique, by = "species")

# guild_fg column has ectomycorrhizal assignment that I'm interested in 
# Look at the different assignments by guild here 
guild_count_species <- otu_with_traits_species %>%
  count(guild_fg) %>%
  filter(n > 1)
# very very few with any usable assignment at this level 


# 2.Try to match OTUs using just the Genus
otu_with_traits_genus <- ITS_otus %>%
  left_join(Traits_unique, by = "Genus")

# will have many-to-many matching because there are a lot of duplicate genera 

# Look at the different assignments by guild here 
guild_count_genus <- otu_with_traits_genus %>%
  count(guild_fg) %>%
  filter(n > 1)

# Gotta be doing it at the genus level! 


# take genus assignments dataframe and check for variation in assignments within a genus 

# sort to only include those with guild_fg assignments 

otu_with_traits_genus$guild_fg <- as.factor(otu_with_traits_genus$guild_fg)

traits_genus_trim <- otu_with_traits_genus %>%
  filter(!is.na(guild_fg))
# this gives 10,896 rows that have genus and guild assignments 


# Now want to sort down to just guilds that are relevant to ectomycorrhizae
# There are 6 classifications that include 'Ectomycorrhizal' in the name 

traits_genus_ecto <- traits_genus_trim %>%
  filter(guild_fg %in% c("Ectomycorrhizal", "Ectomycorrhizal-Endophyte-Ericoid Mycorrhizal-Litter Saprotroph-Orchid Mycorrhizal",
                         "Ectomycorrhizal-Fungal Parasite", "Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph",
                         "Ectomycorrhizal-Undefined Saprotroph", "Ectomycorrhizal-Wood Saprotroph"))
# 5,432 rows, still with some duplication 


# Find currently duplicated OTUs to double check before filtering to just unique 

dupe_otus <- traits_genus_ecto %>%
  add_count(OTU, name = "otu_count") %>%
  filter(otu_count > 1)
# At most an OTU is duplicated 8 times, at least twice (obvs)

# There is some variation for the species coming from the FungalTraits dataset, but there is no variation in 
# the guild_fg assignment at the genus level, so there should be no problem with trimming duplicate OTUs. 

# Keep distinct OTUs for the final dataset 

final_otu_traits <- traits_genus_ecto %>%
  distinct(OTU, .keep_all = TRUE)

# 2,140 retained 

# Check it out 

final_otu_traits$Genus <- as.factor(final_otu_traits$Genus)


# Save dataset for trait assignments from the literature 
write.csv(final_otu_traits, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/EM_FungalTraits_final.csv")


## Conclusion was that while more ASV's were given assignments, 15 fewer genera had assignments than what 
# was done using FunGuild, so the coverage was much poorer. 

# This is bring retained as a comparison, but the original FunGuild assignments were used for all downsteam steps 

