# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Processing raw amplicon reads
# Author: Geoffrey Zahn
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dada2 v 1.32.0
#                     decontam v 1.24.0
#                     phyloseq v 1.48.0
#                     purrr v 1.0.4
#                     Biostrings v 2.72.1
#                     patchwork v 1.3.0
#                     dplyr v 1.1.4
#                     stringr v 1.5.1
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inference, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(dplyr); packageVersion("dplyr")
library(stringr); packageVersion("stringr")

# PARSE FILE PATHS ####

# File parsing - 

path <- "/Users/mckinleynevins/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Data/240223_Nevins_AMF_demultiplxed_fastqs/cutadapt" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
#using the cutadapted AM files 


#these add structures to save the filtered files into the filtN folder 
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Changed the pattern to match the naming convention for my forward and reverse reads 
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_001.fastq.gz")) # make pattern match your FWD reads
rns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R2_001.fastq.gz")) # make pattern match your REV reads


# The following line splits the filename on "_" and keeps the first element
# This turns the same names into just the location-species-replicate, for example (A-ABAM-01)
sample.names <- unlist(map(strsplit(basename(fns), "_"), 1)) # this pulls out just the basename

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files

# can check in on any of the forward reads for any sample 
p1 <- plotQualityProfile(fns[150]) + ggtitle("Example forward reads")
p2 <- plotQualityProfile(rns[150]) + ggtitle("Example reverse reads")

# display and save the plots
p1 / p2
# ggsave("./output/figs/unfiltered_quality_plots.png",dpi=500,height = 6,width = 6)

# FILTER AND TRIM ####

# only using forward reads after this because that's the standard practice with this primer set. 
# See Davison et al. 2012 as a reference for the methodological choice. 

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))

# this is the actual quality control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2, and will be removed
                     maxEE=c(2), # refers to the maximum expected errors allowed for forward and reverse
                     truncQ=2, # special value denoting "end of good quality sequence" (optional), regular # for phred score
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     truncLen = c(200), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=4) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./trackreads_cutadapt_test_AM.RDS")

# Did any samples have NO reads pass the filtration step?

length(fns);length(filts_f) # will be the same length if all samples had some passing reads
# 219 and 219. Looks like some blanks had none pass, but not sure why that isn't reflected here. 

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))

# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[150]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[141])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4


# ggsave("./output/figs/filtered_quality_comparison.png",dpi=500,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
set.seed(123) # "random" seed for reproducibility
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errF,"./errF_cutadapt_test_AM.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./error_model_F_cutadapt_test_AM.png",dpi=500,height = 6,width = 6)

# dereplication - removes duplicate sequences to make the computation faster
derepF <- derepFastq(filts_f, verbose=TRUE)

# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derepF) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts_f), "_filt"), 1))
}

##the next step may be the most time consuming one!! Depends on the length or reads, the diversity of 
# samples, and the processing power available 


# SAMPLE INFERRENCE #### - giving error profiles, having dada correct any sequencing mistakes. Makes assumption that most abundant read is only real one
set.seed(123)
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, MAX_CONSIST = 20, pool = "pseudo") # set multithread = FALSE on Windows

# MAKE SEQUENCE TABLE ####

seqtab <- makeSequenceTable(dadaFs)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 935 bimeras out of 6172 input sequences.
saveRDS(seqtab.nochim,"./seqtab.nochim_cutadapt_test_AM.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./AMF_read_counts_at_each_step_cutadapt.csv", row.names = TRUE)


# IMPORT METADATA ####
meta <- read_csv("./am_community_metadata_clean.csv")
meta <- as.data.frame(meta)
row.names(meta) <- meta$SampleNumber
meta <- as_tibble(meta)

df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                 SampleNumber=row.names(seqtab.nochim))
df2 <- left_join(meta,df,by="SampleNumber")
df2 <- as.data.frame(df2)
row.names(df2) <- df2$SampleNumber
meta <- df2[row.names(seqtab.nochim),]
identical(row.names(meta),row.names(seqtab.nochim)) #TRUE


# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_asvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_asvs]

# remove singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./seqtab.nochim.clean_cutadapt_AM.RDS")

seqtab.nochim <- readRDS("./seqtab.nochim.clean_cutadapt_AM.RDS")

# Find and remove contaminants ####

# this says control, but what it really means is the blanks we ran with each extraction 
# the 'control' that will be left in my data are the actual controls I did extractions 
# for using the MycoBloom 
meta$Control <- meta$Host_ID == "Blank"
meta <- meta[meta$SampleNumber %in% row.names(seqtab.nochim),] 

contams = isContaminant(seqtab.nochim, neg = meta$Control, normalize = TRUE)

table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./likely_contaminants_cutadapt_AM.csv", row.names = TRUE)

# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Control == FALSE,]
meta = meta[meta$Control == FALSE,]


################### SAVE FILES FOR LATER USE ##########

#save modified meta file 
write.csv(meta, file = "./AM_meta_final_cutadapt_AM.csv", row.names = TRUE)

#save seqtab with contaminants removed. Fully ready for taxonomy - 5,235 ASV's
saveRDS(seqtab.nochim,"./AM.seqtab.nochim.final_cutadapt_AM.RDS")

# ASSIGN TAXONOMY #### - do these steps on Kamiak since it's so intensive 

# downloaded the full Maarjam database, relevant for environmental samples of AMF communities 
# using the new version of the database 
# this is stored in the 'data' file, so switch the working directory 
maarjam.ref <- "./Eukaryome_General_SSU_v1.8_reformatted_VTX.fasta.gz"
taxa <- assignTaxonomy(seqtab.nochim, maarjam.ref, multithread=4, tryRC=TRUE, verbose=TRUE)

# Save completed taxonomy file
saveRDS(taxa, file = "./AM_Taxonomy_final_2025.RDS")

##################################################
#load in taxonomy file from Kamiak

AM_taxa <-read.table("AM_Taxonomy_final_2025.txt",header=TRUE,sep=" ")
#assignments made for all 5,235 ASV's

taxa <- as.matrix(AM_taxa)

# inspect taxonomy
AM_taxa.print <- AM_taxa # Removing sequence rownames for display only
rownames(AM_taxa.print) <- NULL
head(AM_taxa.print)

library('phyloseq')

#read in final and clean sequence table for AM
AM_seqtab.nochim <- readRDS("./AM.seqtab.nochim.final_cutadapt_AM.RDS")

#read in modified meta file 
AM_meta <- read.csv("./am_community_metadata_clean.csv")


# Hand off to Phyloseq ####
otu <- otu_table(AM_seqtab.nochim, taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(AM_meta)
row.names(met) <- AM_meta$SampleNumber # more rows here than what gets pulled into the phyloseq 
                                        # object because the controls and blanks are gone 


AM_ps <- phyloseq(otu,met,tax)
saveRDS(AM_ps,"./AM_ps_not-cleaned_2025.RDS")
AM_ps

