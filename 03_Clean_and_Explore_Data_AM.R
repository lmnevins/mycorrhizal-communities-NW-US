# -----------------------------------------------------------------------------#
# Clean AM phyloseq object and prepare for FUNGuild 
# Original Author: Geoffrey Zahn
# Modified by: L. McKinley Nevins 
# March 10, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     ShortRead v 1.62.0
#                     phyloseq v 1.48.0
#                     Biostrings v 2.72.1
#                     microViz v 0.12.4
#                     microbiome v 1.26.0
#                     ComplexHeatmap v 2.20.0
#                     
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#  Remove plant, bacteria, and other non-fungal sequences from the data.        #
#  Prepare the otu table to be processed through the FUNGuild database to make  #
#  guild assignments.                                                           #
#                                                                               #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####
BiocManager::install(c("microbiome", "ComplexHeatmap"), update = FALSE)

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(microViz); packageVersion("microViz")
library(microbiome); packageVersion("microbiome")
library(ComplexHeatmap); packageVersion("ComplexHeatmap")

# load phyloseq object 

ps_AM <- readRDS("./AM_ps_not-cleaned_2025.RDS") # 

# Find non-fungi ####
ps_nonfungi <- subset_taxa(ps_AM, Kingdom != "k__Fungi")

#look at non-fungi 
non_fung <- tax_table(ps_nonfungi) %>% as("matrix") %>% as.data.frame()
#4,085 sequences that are plants, insects, etc. 

#Check what was in the controls 
ps_AM_control <- ps_filter(ps_AM, Host_ID == "Control")

#look at controls 
controls <- tax_table(ps_AM_control) %>% as("matrix") %>% as.data.frame()

#filter to remove controls - using the microViz package
ps_AM <- ps_filter(ps_AM, Host_ID != "Control")

check <- otu_table(ps_AM) %>% as("matrix") %>% as.data.frame()

count <- rowSums(check) %>% as.data.frame()

#filter to remove one sample that didn't have any fungi identified - South CONU 01
# waiting to do this, there may be others after the cleaning for non-fungal taxa 
# that also need to be trimmed. 
#ps_AM <- ps_filter(ps_AM, SampleNumber != "S-CONU-01_FWD")


# quick plot to look at kingdom-level taxonomy
ps_AM %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")

#ggsave("./Phylogeny_Outputs/Figs/AM_Kingdom-Level_Taxonomic_Proportions.png",dpi=300) # save figure for later use

# samples are mostly fungi, but have a lot of viridiplantae and mitochondria mixed in 

# same plot, but non-fungi, for sanity check
ps_nonfungi %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")

# REMOVE NON-Fungi, CHLOROPLASTS, MITOCHONDRIA, and empty samples/taxa ####
ps_AM <- subset_taxa(ps_AM, Kingdom == "k__Fungi")
ps_AM_tax <- tax_table(ps_AM) %>% as("matrix") %>% as.data.frame()

# 629 ASV's 

############
# This is something that could be returned to if I want to 
# exclude more samples or taxa that don't have many reads, etc. 
# If they don't seem representative of the real community 
############


#check to make sure samples still have reads 
check <- otu_table(ps_AM) %>% as("matrix") %>% as.data.frame()

count <- rowSums(check) %>% as.data.frame()

# look at which taxa still have reads 
count2 <- colSums(check) %>% as.data.frame()

# a portion of taxa are quite rare 

ps_AM <- subset_taxa(ps_AM, taxa_sums(ps_AM) > 0)
ps_AM <- subset_samples(ps_AM, sample_sums(ps_AM) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps_AM)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./Phylogeny_Outputs/AM_OTU_reference_sequences_2025.RDS")

# Save RDS object for cleaned up Phyloseq object
saveRDS(ps_AM, file = "./Phylogeny_Outputs/AM_clean_phyloseq_2025.RDS")
ps_AM


##########################################################################################
# save OTU table separately to be processed by FUNGuild

ps_AM_otu <- otu_table(ps_AM) %>% as("matrix") %>% as.data.frame()

# for the FUNGuild format, the samples need to be in the header and the OTU's need 
# to be individual rows 

#transpose the columns and rows of the otu table 
ps_AM_otu <- t(ps_AM_otu)
#worked


# then there needs to be a column at the end called 'taxonomy', which has the data assigned 
# from Maarjam or another database 

# grab taxa assignments 
ps_AM_tax <- tax_table(ps_AM) %>% as("matrix") %>% as.data.frame()


# now just need to merge the taxonomic ranks into one column of ; separated values 
ps_AM_tax$taxonomy <- paste(ps_AM_tax$Kingdom, ps_AM_tax$Phylum, ps_AM_tax$Class, ps_AM_tax$Order,
                             ps_AM_tax$Family, ps_AM_tax$Genus, ps_AM_tax$Species, sep=";")

#subset to get just the otu and taxonomy columns (OTU's are rownames right now)
ps_AM_tax <- select(ps_AM_tax, taxonomy)

# join the OTU object and the taxonomy object by columns that match up OTU rows
AM_table <- cbind(ps_AM_otu, ps_AM_tax)

#get rid of rownames 
AM_table <- tibble::rownames_to_column(AM_table, "OTU")

# dataframe is now organized to have an OTU column leading, followed by one column for 
# each sample, and a finally taxonomy column with each rank separated by ';' 

#save dataframe as a txt file 
write.csv(AM_table, "./FUNGuild/AM_OTU_table_2025.csv")

#saves with row numbers in the first column. After it saves, open it up in Excel just 
#to double check it, and delete that first column. 




