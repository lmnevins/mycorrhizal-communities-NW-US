# -----------------------------------------------------------------------------#
# Clean ITS phyloseq object and prepare for FUNGuild 
# Original Author: Geoffrey Zahn
# Modified by: L. McKinley Nevins 
# July 5, 2024
# Software versions:  R v 4.2.1
#                     tidyverse v 1.3.2
#                     ShortRead v 
#                     phyloseq v 1.42.0
#                     Biostrings v 2.66.0

# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#  Remove plant, bacteria, and other non-fungal sequences from the data.        #
#  Prepare the otu table to be processed through the FUNGuild database to make  #
#  guild assignments.                                                           #
#                                                                               #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(microViz)

# load phyloseq object with phylogenetic tree ####

ps_ITS <- readRDS("./Phylogeny_Outputs/ITS_ps_not-cleaned_w_tree.RDS") # change to non-phylogeny stuff

# Find non-fungi ####
ps_nonfungi <- subset_taxa(ps_ITS, Kingdom != "k__Fungi")

#look at non-fungi 
non_fung <- tax_table(ps_nonfungi) %>% as("matrix") %>% as.data.frame()

#Check what was in the controls 
ps_ITS_control <- ps_filter(ps_ITS, Host_ID == "Control")

#filter to remove controls - using the microViz package
ps_ITS <- ps_filter(ps_ITS, Host_ID != "Control")

#filter to remove non-fungi from samples 
ps_ITS_control <- subset_taxa(ps_ITS_control, Kingdom == "k__Fungi")

ps_ITS_control %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")

ps_ITS_control %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Phylum")

ps_ITS_control %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Genus")


# quick plot to look at kingdom-level taxonomy
ps_ITS %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")

ggsave("./Phylogeny_Outputs/Figs/ITS_Kingdom-Level_Taxonomic_Proportions.png",dpi=300) # save figure for later use

# same plot, but non-fungi, for sanity check
ps_nonfungi %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Phylum")

# REMOVE NON-FUNGAL - Metazoa, Rhizaria, Rhodoplantae, Stramenopila, 
# Viridiplantae - and empty samples/taxa ####
ps_ITS <- subset_taxa(ps_ITS, Kingdom == "k__Fungi")
ps_ITS_tax <- tax_table(ps_ITS) %>% as("matrix") %>% as.data.frame()

#I don't actually think I need to do this, because the previous subset puts only Fungi 
#into the dataset 
#ps_ITS <- subset_taxa(ps_ITS,Class != "Chloroplast")

ps_ITS <- subset_taxa(ps_ITS, taxa_sums(ps_ITS) > 0)
ps_ITS <- subset_samples(ps_ITS, sample_sums(ps_ITS) > 0)


# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps_ITS)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./Phylogeny_Outputs/ITS_ASV_reference_sequences.RDS")

# Save RDS object for cleaned up Phyloseq object
saveRDS(ps_ITS, file = "./Phylogeny_Outputs/ITS_clean_phyloseq_object.RDS")
ps_ITS

##########################################################################################
# save OTU table separately to be processed by FUNGuild
 
ps_ITS_otu <- otu_table(ps_ITS) %>% as("matrix") %>% as.data.frame()

# for the FUNGuild format, the samples need to be in the header and the OTU's need 
# to be individual rows 

#transpose the columns and rows of the otu table 
ps_ITS_otu <- t(ps_ITS_otu)
#worked


# then there needs to be a column at the end called 'taxonomy', which has the data assigned 
# from UNITE or another database 

# grab taxa assignments 
ps_ITS_tax <- tax_table(ps_ITS) %>% as("matrix") %>% as.data.frame()


# now just need to merge the taxonomic ranks into one column of ; separated values 
ps_ITS_tax$taxonomy <- paste(ps_ITS_tax$Kingdom, ps_ITS_tax$Phylum, ps_ITS_tax$Class, ps_ITS_tax$Order,
                       ps_ITS_tax$Family, ps_ITS_tax$Genus, ps_ITS_tax$Species, sep=";")

#subset to get just the otu and taxonomy columns (OTU's are rownames right now)
ps_ITS_tax <- select(ps_ITS_tax, taxonomy)

# join the OTU object and the taxonomy object by columns that match up OTU rows
ITS_table <- cbind(ps_ITS_otu, ps_ITS_tax)

#get rid of rownames 
ITS_table <- tibble::rownames_to_column(ITS_table, "OTU")

# dataframe is now organized to have an OTU column leading, followed by one column for 
# each sample, and a finally taxonomy column with each rank separated by ';' 

#save dataframe as a txt file 
write.csv(ITS_table, "./FUNGuild/EM_OTU_table.csv")

#saves with rownumbers in the first column. After it saves, open it up in Excel just 
#to doublecheck it, and delete that first column. 





