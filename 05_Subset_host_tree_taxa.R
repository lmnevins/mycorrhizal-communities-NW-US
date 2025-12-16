# -----------------------------------------------------------------------------#
# Take EM and AM communities from FungalTraits and FunGuild processing to 
# guilds and subset to target host tree taxa
# Original Author: L. McKinley Nevins 
# March 11, 2025 for AM
# December 15, 2025 for EM
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
#  Process the FUNGuild communities to target just host tree taxa of interest.  #
#  CONU and PIMO do not have full coverage across the study sites and need to   #
#  be removed before all of the community analyses. This will leave in the AM   #
#  community: ALRU, TABR, AND THPL. In the EM community: ABAM, ABGR, ABPR,      #
#  ALRU, PSME, TABR, TSHE. The resulting files should be used for all of the    #
#  downstream analyses.                                                         #
#                                                                               #
#################################################################################

################ -- 
# (1) DATA PREP
################ -- 

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/"
setwd(wd)


##phyloseq objects:
# Load phyloseq object produced from filtered FUNGuild guilds to target just AM fungi ####
ps_AM <- readRDS("./Phylogeny_Outputs/AM_funguild_final_2025.RDS")

# Load phyloseq object produced from filtered FUNGuild guilds to target just EM fungi ####
ps_EM <- readRDS("./Phylogeny_Outputs/EM_funguild_final_2025.RDS")


##environmental data:
# Load in the sample data file containing all of the environmental data for the AM community 
# of just ALRU, TABR, and THPL that were retained with reads 
# these are in the 'FINAL' file
env_AM <- read.csv(file = "./FINAL/AM_subset_sample_data_2025.csv", row.names = 1)

# Load in the sample data file containing all of the environmental data for the EM community 
# for the subset of 7 host tree species, excluding THPL as it is only an AM host  

# This results in 372 host trees 
env_EM <- read.csv(file = "./FINAL/EM_subset_sample_data_2025.csv", row.names = 1)



###Subsetting the host trees to get the final communities for all downstream analyses 
###################
#AM Community 

#preparation of the tax_table from the full AM phyloseq object
tax_AM <- tax_table(ps_AM)

tax_AM <- as.data.frame(tax_AM)

#convert to phyloseq compatible object 
tax_AM <- phyloseq::tax_table(as.matrix(tax_AM))

otu_AM <- otu_table(ps_AM)

otu_AM <- as.data.frame(otu_AM)

write.csv(otu_AM, "./FINAL/AM_subset_otu_2025.csv")

#read in the otu file of just the subset of AM target trees moving forward - ALRU, TABR, THPL
#this was generated from the otu table above, excluding CONU
AM_subset_otu <- read.csv("./FINAL/AM_subset_otu_2025.csv")

#need the samples to be row names
AM_subset_otu <- data.frame(AM_subset_otu[,-1], row.names=AM_subset_otu[,1])

#convert to phyloseq compatible object 
AM_subset_otu <- phyloseq::otu_table(as.matrix(AM_subset_otu), taxa_are_rows = F)


###create subset phyloseq object 
ps_AM_final <- phyloseq(otu_table(AM_subset_otu), tax_table(tax_AM), sample_data(env_AM))

#inspect
# number of taxa - 438
ntaxa(ps_AM_final)

# number of samples - 130
nsamples(ps_AM_final)

asv <- otu_table(ps_AM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

sample_count <- rowSums(asv)

sample_count <- as.data.frame(sample_count)

# All trees have at least one ASV

asv_count <- colSums(asv)

asv_count <- as.data.frame(asv_count) #some ASV's no longer present in the samples 
# these were essentially the ones that were just associated with CONU

#remove taxa that aren't present in this new subset 
ps_AM_final <- subset_taxa(ps_AM_final, taxa_sums(ps_AM_final) > 0)

# number of taxa - now 399
ntaxa(ps_AM_final)

#SAVE THE FINAL PHYLOSEQ OBJECT 
saveRDS(ps_AM_final, file = "./FINAL/AM_phyloseq_final_2025.RDS")

#now I can save the tax table and look at the community composition 
AM_tax <- tax_table(ps_AM_final)

AM_tax <- as.data.frame(AM_tax)

write.csv(AM_tax, "./FINAL/AM_subset_tax_2025.csv")

############################################################################
#EM Community 

#preparation of the tax_table from the full EM phyloseq object
tax_EM <- tax_table(ps_EM)

tax_EM <- as.data.frame(tax_EM)

#convert to phyloseq compatible object 
tax_EM <- phyloseq::tax_table(as.matrix(tax_EM))


#read in the otu file of just the subset of EM target trees moving forward - 
#7 of the 10 species excluding CONU, PIMO, and THPL
#this was generated from the otu table from the EM funguild phyloseq object 
EM_subset_otu <- read.csv("./FINAL/EM_subset_otu_2025.csv")

#need the samples to be row names
EM_subset_otu <- data.frame(EM_subset_otu[,-1], row.names=EM_subset_otu[,1])

#convert to phyloseq compatible object 
EM_subset_otu <- phyloseq::otu_table(as.matrix(EM_subset_otu), taxa_are_rows = F)


###create subset phyloseq object 
ps_EM_final <- phyloseq(otu_table(EM_subset_otu), tax_table(tax_EM), sample_data(env_EM))


#inspect
# number of taxa - 2,343
ntaxa(ps_EM_final)

# number of samples - 372
nsamples(ps_EM_final)

asv <- otu_table(ps_EM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

sample_count <- rowSums(asv) %>% as.data.frame()
#all trees with at least some ASV's

asv_count <- colSums(asv) %>% as.data.frame()
#some ASV's no longer present in the samples 

#remove taxa that aren't present in this new subset 
ps_EM_final <- subset_taxa(ps_EM_final, taxa_sums(ps_EM_final) > 0)

# number of taxa - now 1,870
ntaxa(ps_EM_final)

#SAVE THE FINAL PHYLOSEQ OBJECT 
saveRDS(ps_EM_final, file = "./FINAL/EM_phyloseq_final_2025.RDS")


#now I can save the tax table and look at the community composition 
EM_tax <- tax_table(ps_EM_final)

EM_tax <- as.data.frame(EM_tax)

write.csv(EM_tax, "./FINAL/EM_subset_tax_2025.csv")

