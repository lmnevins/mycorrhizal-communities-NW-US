# -----------------------------------------------------------------------------#
# Process the assigned FUNGuild objects for EM and AM communities
# Original Author: L. McKinley Nevins 
# July 8, 2024 for ITS (EM)
# March 11, 2025 for AM
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


# Load FUNGuild output files ####

AM_funguild <- read.csv("./FUNGuild/AM_OTU_table_2025.guilds_matched.csv")
# 558 OTU's with assignments 


ITS_funguild <- read.csv("./FUNGuild/ITS_otu_table.guilds_matched_sorted.csv")
# 11,984 OTU's with assignments 


### Make guild into factor for both so it can be sorted ####

AM_funguild$Guild <- as.factor(AM_funguild$Guild)

#make list to inspect 
AM_guilds <- as.data.frame(AM_funguild$Guild)
# many guilds 


ITS_funguild$Guild <- as.factor(ITS_funguild$Guild)

#make list to inspect 
ITS_guilds <- as.data.frame(ITS_funguild$Guild)
# 171 guild categories 

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

# 27 guilds that are ectomycorrhizal 
EM_funguild <- ITS_funguild [ITS_funguild$Guild %in% c("|Ectomycorrhizal|", "|Ectomycorrhizal|-Fungal Parasite",
                                                        "|Ectomycorrhizal|-Orchid Mycorrhizal-Root Associated Biotroph",
                                                        "|Ectomycorrhizal|-Plant Saprotroph", "|Ectomycorrhizal|-Undefined Saprotroph",
                                                        "Ectomycorrhizal", "Ectomycorrhizal-|Endophyte|-Ericoid Mycorrhizal-Orchid Mycorrhizal-Plant Saprotroph-Undefined Saprotroph",
                                                        "Ectomycorrhizal-|Endophyte|-Ericoid Mycorrhizal-Orchid Mycorrhizal-Undefined Saprotroph",
                                                        "Ectomycorrhizal-|Plant Saprotroph|", "Ectomycorrhizal-|Plant Saprotroph|-Undefined Saprotroph-Wood Saprotroph",
                                                        "Ectomycorrhizal-|Plant Saprotroph|-Wood Saprotroph", "Ectomycorrhizal-|Undefined Saprotroph|-Wood Saprotroph",
                                                        "Ectomycorrhizal-Endomycorrhizal-Orchid Mycorrhizal-|Plant Pathogen|-Plant Saprotroph-Undefined Saprotroph",
                                                        "Ectomycorrhizal-Endophyte-|Undefined Saprotroph|", "Ectomycorrhizal-Endophyte-Plant Pathogen-|Plant Saprotroph|",
                                                        "Ectomycorrhizal-Endophyte-Plant Pathogen-|Plant Saprotroph|-Wood Saprotroph", 
                                                        "Ectomycorrhizal-Ericoid Mycorrhizal-Plant Pathogen-|Plant Saprotroph|-Undefined Saprotroph-Wood Saprotroph", 
                                                        "Ectomycorrhizal-Fungal Parasite", "Ectomycorrhizal-Fungal Parasite-|Undefined Saprotroph|-Undefined Symbiotroph",
                                                        "Ectomycorrhizal-Fungal Parasite-Plant Pathogen-Wood Saprotroph", "Ectomycorrhizal-Fungal Parasite-Plant Saprotroph-Wood Saprotroph",
                                                        "Ectomycorrhizal-Fungal Parasite-Soil Saprotroph-Undefined Saprotroph", "Ectomycorrhizal-Fungal Pathogen-Undefined Saprotroph", 
                                                        "Ectomycorrhizal-Lichen Parasite-Lichenized-Plant Pathogen", "Ectomycorrhizal-Undefined Saprotroph", 
                                                        "Ectomycorrhizal-Undefined Saprotroph-Wood Saprotroph", "Endomycorrhizal-Plant Pathogen-Undefined Saprotroph"), ]

# 2,542 OTU's remaining after filtering 
# save the dataframe 
write.csv(EM_funguild, "./FUNGuild/EM_guilds_final.csv", row.names = FALSE)


### Load the funguild tables back into Phyloseq objects 

# need to have an OTU table and a taxa table to load back into the phyloseq object with the sample data 

#OTU table has each OTU as a row and the samples as columns, need to transpose and have the samples as row names

#Taxa table has each OTU as a row, and each of the taxonomic ranks as a separate column. Could separate 
#out the ranks that I merged together into the taxonomy column, or I could merge the OTU's and it should keep 
#only the taxonomy for the OTUs present in the funguild dataset 

# made the tax table outside of R from the csv file 
# made the otu table

###AM####

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

###EM####

#tax table 
EM_funguild_tax <- read.csv("./FUNGuild/EM_funguild_tax.csv")

#need the OTUs to be row names with no OTU column 
EM_funguild_tax <- data.frame(EM_funguild_tax[,-1], row.names=EM_funguild_tax[,1])
#DONE

####
#otu table 
EM_funguild_otu <- read.csv("./FUNGuild/EM_funguild_otu.csv")

#need the samples to be row names
EM_funguild_otu <- data.frame(EM_funguild_otu[,-1], row.names=EM_funguild_otu[,1])
#DONE

#sample data 
#shouldn't need any adjustment right now because it's all the same samples as previously
EM_sample_data <- (sample_data(ps_ITS))
#DONE

#put it all together 

EM_funguild_otu <- phyloseq::otu_table(as.matrix(EM_funguild_otu), taxa_are_rows = F)

EM_funguild_tax <- phyloseq::tax_table(as.matrix(EM_funguild_tax))

#sample data is already a phyloseq sample data object 

ps_EM_final <- phyloseq(otu_table(EM_funguild_otu), tax_table(EM_funguild_tax), sample_data(EM_sample_data))

#can use these to check that the names are reading in the same, if an error message arises 
taxa_names(EM_funguild_otu)

taxa_names(EM_funguild_tax)

#filter to remove one sample that didn't have any fungi identified - North ABPR 09
ps_EM_final <- ps_filter(ps_EM_final, SampleNumber != "N-ABPR-09_FWD")

##Saving the final object
saveRDS(ps_EM_final, file = "./Phylogeny_Outputs/EM_funguild_final.RDS")



