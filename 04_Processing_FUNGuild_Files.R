# -----------------------------------------------------------------------------#
# Process the assigned FUNGuild objects for EM and AM communities
# Original Author: L. McKinley Nevins 
# March 11, 2025 for AM
# December 14, 2025 for ITS (EM)
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     phyloseq v 1.48.0
#                     vegan v 2.6.10
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")


#################################################################################
#                               Main workflow                                   #
#    Process the FUNGuild outputs to make sure they are sorted properly for     #
#    later analyses. Removing non-target taxa.                                  #
#                                                                               #
#################################################################################

################ -- 
# (1) DATA PREP
################ -- 

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/"
setwd(wd)


## FunGuild was run through the command line for the AM and ITS data. The AM community 
# was everything filtered as AM ASVs. The ITS data was everything that had been 
# previously run through FungalTraits with no successful guild assignment, but which 
# had genus level info. 

# Load FUNGuild output files ####

AM_funguild <- read.csv("./FUNGuild/AM_OTU_table_2025.guilds_matched.csv")
# 558 ASVs with assignments, out of 629 provided to FUNGuild

# From the original FUNGuild output file, the original txt file is opened and 
# saved as a csv
ITS_funguild <- read.csv("./FUNGuild/ITS_OTU_table_2025.guilds_matched.csv")
# 7,565 ASVs with assignments, out of the 13,458 provided to FUNGuild, of which 
# 11,132 had matching taxonomy records in the database 


### Make guild into factor for both so it can be sorted ####

AM_funguild$Guild <- as.factor(AM_funguild$Guild)

#make list to inspect 
AM_guilds <- as.data.frame(AM_funguild$Guild)
# many guilds 


ITS_funguild$Guild <- as.factor(ITS_funguild$Guild)

#make list to inspect 
ITS_guilds <- as.data.frame(ITS_funguild$Guild)
# 138 guild categories 


#################################################################################

###################### -- 
# (2) FILTER GUILDS 
###################### -- 

### Filtering to remove any non-AM and EM OTU's ######

##AM#####

# variety of non-AM guilds, so it'll be easier to select the ones I want and to 
# exclude everything else 

# "|Arbuscular Mycorrhizal|-Endophyte", "Arbuscular Mycorrhizal"

#selecting all guilds that are these two classifications 
AM_funguild <- AM_funguild [AM_funguild$Guild %in% c("|Arbuscular Mycorrhizal|-Endophyte", "Arbuscular Mycorrhizal"), ]

# 438 remaining 
# excluded 120 that were non AMF

# save the dataframe 
write.csv(AM_funguild, "./FUNGuild/AM_guilds_final_2025.csv", row.names = FALSE)

##EM########

# There are a ton of non-ectomycorrhizal guilds, so it will be faster to select the guilds that 
# are EM

# 20 guilds that are ectomycorrhizal 
EM_funguild <- ITS_funguild [ITS_funguild$Guild %in% c("|Ectomycorrhizal|","|Ectomycorrhizal|-Orchid Mycorrhizal-Root Associated Biotroph",
                                                        "|Ectomycorrhizal|-Plant Saprotroph", "|Ectomycorrhizal|-Undefined Saprotroph","Ectomycorrhizal",
                                                        "Ectomycorrhizal-|Plant Saprotroph|", "Ectomycorrhizal-|Plant Saprotroph|-Undefined Saprotroph-Wood Saprotroph",
                                                        "Ectomycorrhizal-|Plant Saprotroph|-Wood Saprotroph",
                                                        "Ectomycorrhizal-Endomycorrhizal-Orchid Mycorrhizal-|Plant Pathogen|-Plant Saprotroph-Undefined Saprotroph",
                                                        "Ectomycorrhizal-Endophyte-Plant Pathogen-|Plant Saprotroph|-Wood Saprotroph", 
                                                        "Ectomycorrhizal-Ericoid Mycorrhizal-Plant Pathogen-|Plant Saprotroph|-Undefined Saprotroph-Wood Saprotroph", 
                                                        "Ectomycorrhizal-Fungal Parasite",
                                                        "Ectomycorrhizal-Fungal Parasite-Plant Pathogen-Wood Saprotroph", "Ectomycorrhizal-Fungal Parasite-Plant Saprotroph-Wood Saprotroph",
                                                        "Ectomycorrhizal-Fungal Parasite-Soil Saprotroph-Undefined Saprotroph", "Ectomycorrhizal-Fungal Pathogen-Undefined Saprotroph", 
                                                        "Ectomycorrhizal-Lichen Parasite-Lichenized-Plant Pathogen", "Ectomycorrhizal-Undefined Saprotroph", 
                                                        "Ectomycorrhizal-Undefined Saprotroph-Wood Saprotroph", "Endomycorrhizal-Plant Pathogen-Undefined Saprotroph"), ]

# 422 OTU's remaining after filtering 


# Check the confidence ratings 
EM_funguild$Confidence_Ranking <- as.factor(EM_funguild$Confidence.Ranking)

summary(EM_funguild$Confidence_Ranking)

# 43 highly probable, 247 probable, 132 possible 

# Filter to remove any with 'Possible' assignments 
EM_funguild <- EM_funguild %>%
  filter(Confidence.Ranking != "Possible")

# 290 remaining ASVs


# Split taxonomy column to get separate ranks 
EM_funguild <- EM_funguild %>%
  separate(
    col = taxonomy,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    fill = "right",   
    remove = FALSE    
  )


# Fix name format for genus and species 

EM_funguild$Genus <- sub("^g__", "", EM_funguild$Genus)

EM_funguild$Species <- sub("^s__", "", EM_funguild$Species)


# save the dataframe 
write.csv(EM_funguild, "./FUNGuild/EM_guilds_final_2025.csv", row.names = FALSE)


# Filter columns to be able to merge with fungaltraits data down below 
EM_funguild <- dplyr::select(EM_funguild, -c(Taxon, Taxon.Level, Trophic.Mode, Growth.Morphology, 
                                             Notes, Trait, Citation.Source, Confidence.Ranking))


#################################################################################

######################## -- 
# (3) MERGE ASSIGNMENTS  
######################## -- 


# Merge the 290 ASVs with assignments here with the 2,053 ASVs with assignments from 
# FungalTraits to get the final 


# Load in ASVs with assignments from FungalTraits 

EM_fungaltraits <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/EM_FungalTraits_final.csv", 
                            row.names = 1)


# Put a few pieces of data into columns with compatible names 
EM_fungaltraits$Guild <- EM_fungaltraits$guild_fg

EM_fungaltraits$Confidence_Ranking <- EM_fungaltraits$confidence_fg

# Can trim extra columns we don't need anymore 
EM_fungaltraits <- dplyr::select(EM_fungaltraits, -c(species.x, species.y, speciesMatched,
                                                   ifungorum_number, confidence_fg, growth_form_fg, guild_fg,
                                                   higher_clade, substrate, taxonomic_level_fg, trait_fg, 
                                                   trophic_mode_fg))

# This just leaves the combined and separate taxonomy columns, Guild, and Confidence.Ranking 


# Fix column name format 
names(EM_funguild) <- gsub(pattern = "\\.", replacement = "-", x = names(EM_funguild))

names(EM_fungaltraits) <- gsub(pattern = "\\.", replacement = "-", x = names(EM_fungaltraits))


## Check that the columns in the funguild and fungaltraits dataframes are in the same order 
colnames(EM_funguild)
colnames(EM_fungaltraits)

colnames(EM_funguild) == colnames(EM_fungaltraits)

# Names and order match, so can stack 


# Stack dataframes 
EM_combined <- bind_rows(EM_funguild, EM_fungaltraits)


# Check if there are any duplicate ASVs, and remove them 
EM_combined <- EM_combined %>%
  distinct(OTU, .keep_all = TRUE)

# No change, which is to be expected based on all of the previous filtering we did. 


#################################################################################

############################### -- 
# (4) CREATE PHYLOSEQ OBJECTS  
############################### -- 


### Load the final funguild tables back into Phyloseq objects 

# need to have an OTU table and a taxa table to load back into the phyloseq object with the sample data 

#OTU table has each OTU as a row and the samples as columns, need to transpose and have the samples as row names

#Taxa table has each OTU as a row, and each of the taxonomic ranks as a separate column. Could separate 
#out the ranks that I merged together into the taxonomy column, or I could merge the OTU's and it should keep 
#only the taxonomy for the OTUs present in the funguild dataset 

# made the tax table outside of R from the csv file 
# made the otu table

### AM FUNGI 

#tax table 
AM_funguild_tax <- read.csv("./FUNGuild/AM_funguild_tax_2025.csv")

#need the OTUs to be row names with no OTU column 
AM_funguild_tax <- data.frame(AM_funguild_tax[,-1], row.names=AM_funguild_tax[,1])
#DONE

####
#otu table 
AM_funguild_otu <- read.csv("./FUNGuild/AM_funguild_otu_2025.csv", stringsAsFactors = FALSE)

# need to transpose this 
AM_funguild_otu = setNames(data.frame(t(AM_funguild_otu[,-1])), AM_funguild_otu[,1])

# the dashes in the row names have been replaced with periods, which won't match up correctly 

# can grab the row names from the sample data file and load them in as the row names here 

names <- row.names(AM_sample_data)

row.names(AM_funguild_otu) <- names


#DONE

#sample data 
#shouldn't need any adjustment right now because it's all the same samples as previously
AM_sample_data <- (sample_data(ps_AM))
#DONE

#put it all together 

AM_funguild_otu <- phyloseq::otu_table(as.matrix(AM_funguild_otu), taxa_are_rows = F)

AM_funguild_tax <- phyloseq::tax_table(as.matrix(AM_funguild_tax))

#sample data is already a phyloseq sample data object 

ps_AM_final <- phyloseq(otu_table(AM_funguild_otu), tax_table(AM_funguild_tax), sample_data(AM_sample_data))

#can use these to check that the names are reading in the same, if an error message arises 
taxa_names(AM_funguild_otu)

taxa_names(AM_funguild_tax)


# check that samples still have reads 
check <- otu_table(ps_AM_final) %>% as("matrix") %>% as.data.frame()

count <- rowSums(check) %>% as.data.frame()

# a lot of trees now no longer have reads, or an extremely small number 

# look at which taxa still have reads 
count2 <- colSums(check) %>% as.data.frame()
# they all have some, some are just very rare 

#subset the sample trees that have reads (can be harsher later if necessary)
ps_AM_final <- subset_samples(ps_AM_final, sample_sums(ps_AM_final) > 0)

# this cuts it down to 167 trees 

##Saving the final object
saveRDS(ps_AM_final, file = "./Phylogeny_Outputs/AM_funguild_final_2025.RDS")



###EM FUNGI 

#tax table 
# Just want columns for OTU and the taxonomic ranks 
EM_combined_tax <- dplyr::select(EM_combined, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)

write.csv(EM_combined_tax, "./FUNGuild/EM_combined_tax_2025.csv")


#need the OTUs to be row names with no OTU column 
EM_combined_tax <- data.frame(EM_combined_tax[,-1], row.names=EM_combined_tax[,1])

#DONE

####
#otu table 
EM_combined_otu <- dplyr::select(EM_combined, -c(taxonomy, Kingdom, Phylum, Class, Order, Family, Genus, Species,
                                                 Guild, Confidence_Ranking))

#set OTUs as rownames
EM_combined_otu <- data.frame(EM_combined_otu[,-1], row.names=EM_combined_otu[,1])

# transpose so samples are the rows 
EM_combined_otu <- t(EM_combined_otu) %>% as.matrix() %>% as.data.frame()


# Fix rownames now because the format change didn't stick 
rownames(EM_combined_otu) <- gsub(pattern = "\\.", replacement = "-", x = rownames(EM_combined_otu))

write.csv(EM_combined_otu, "./FUNGuild/EM_combined_otu_2025.csv")

#DONE

#sample data 

# load in cleaned ITS phyloseq object 
ps_ITS <- readRDS("./Phylogeny_Outputs/ITS_clean_phyloseq_object.RDS")

# This still contains the sample data for all of the ASVs, not just the filtered ones 

# Pull sample data for merging 
EM_sample_data <- (sample_data(ps_ITS))

# Sample data has different name formats, because the otu table has periods in between the 
# first parts of the sample code 

#DONE

#put it all together 


# Make the otu and tax tables above into phyloseq compatible objects 
EM_combined_otu <- phyloseq::otu_table(as.matrix(EM_combined_otu), taxa_are_rows = F)

EM_combined_tax <- phyloseq::tax_table(as.matrix(EM_combined_tax))

#sample data is already a phyloseq sample data object 

ps_EM_final <- phyloseq(otu_table(EM_combined_otu), tax_table(EM_combined_tax), sample_data(EM_sample_data))

#can use these to check that the names are reading in the same, if an error message arises 
taxa_names(EM_combined_otu)

taxa_names(EM_combined_tax)


# Pick up here: check that every ASV has a count, and that every sample has ASVs.


# how many sequences observed in each sample?
seq_counts <- otu_table(ps_EM_final) %>% rowSums() %>% as.data.frame()

# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_EM_final) %>% colSums() %>% as.data.frame()

# All ASVs are present at least once 

#filter to remove one sample that didn't have any fungi identified - North ABPR 09
ps_EM_final <- ps_filter(ps_EM_final, SampleNumber != "N-ABPR-09_FWD")
# 493 samples 


##Saving the final object
saveRDS(ps_EM_final, file = "./Phylogeny_Outputs/EM_funguild_final_2025.RDS")

## -- END -- ##


