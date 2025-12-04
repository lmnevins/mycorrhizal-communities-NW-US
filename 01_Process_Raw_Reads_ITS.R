# -----------------------------------------------------------------------------#
# Processing raw ITS amplicon reads
# Original Author: Geoffrey Zahn
# Modified by: L. McKinley Nevins
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.2
#                     dada2 v 1.18.0
#                     decontam v 1.18.0
#                     phyloseq v 1.42.0
#                     purrr v 1.0.2
#                     Biostrings v 2.66.0
#                     patchwork v 1.1.2.9000
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inference, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")


# PARSE FILE PATHS ####

# File parsing - 

path <- "/Users/mckinleynevins/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Data/240223_Nevins__ITS_demultiplxed/filtN" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
#using the filtN files that removed N bases 


#these add structures to save the filtered files into - this will save inside of the cutadapt file so they are nested
#and it's clear what the order is 
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

#for ITS, quality is much higher on the forward reads 
p1 <- plotQualityProfile(fns[276]) + ggtitle("Example forward reads")
p2 <- plotQualityProfile(rns[276]) + ggtitle("Example reverse reads")

# display and save the plots
p1 / p2
# ggsave("./output/figs/unfiltered_quality_plots.png",dpi=500,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path(path, "filtered", paste0(sample.names, "_REV_filt.fastq.gz"))

# this is the actual quality control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, rns, filts_r, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2, and will be removed
                     maxEE=c(2,2), # refers to the maximum expected errors allowed for forward and reverse
                     truncQ=2, # special value denoting "end of good quality sequence" (optional), regular # for phred score
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     truncLen = c(250,200), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=4) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./trackreads.RDS")

# Did any samples have NO reads pass the filtration step?

length(fns);length(filts_f) # will be the same length if all samples had some passing reads
length(rns);length(filts_r) # will be the same length if all samples had some passing reads

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))

# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[100]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[100])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4


# ggsave("./output/figs/filtered_quality_comparison.png",dpi=500,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
set.seed(123) # "random" seed for reproducibility
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errF,"./errF.RDS")
set.seed(123)
errR <- learnErrors(filts_r, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errR,"./errR.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./error_model_F.png",dpi=500,height = 6,width = 6)
plotErrors(errR, nominalQ=FALSE)
ggsave("./error_model_R.png",dpi=500,height = 6,width = 6)

# dereplication - removes duplicate sequences to make the computation faster
derepF <- derepFastq(filts_f, verbose=TRUE)
derepR <- derepFastq(filts_r, verbose=TRUE)

# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derepF) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts_f), "_filt"), 1))
}


if(identical(unlist(map(strsplit(basename(filts_f), "FWD_filt"), 1)),unlist(map(strsplit(basename(filts_r), "REV_filt"), 1)))){
  names(derepF) <- sample.names
  names(derepR) <- sample.names
} else {
  stop("Make sure fwd and rev files are in same order!")
}  


##the next step may be the most time consuming one!! Depends on the length or reads, the diversity of 
# samples, and the processing power available 


# SAMPLE INFERRENCE #### - giving error profiles, having dada correct any sequencing mistakes. Makes assumption that most abundant read is only real one
set.seed(123)
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
set.seed(123)
dadaRs <- dada(derepR, err=errR, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows

# MERGE FWD and REV READS ####
mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=TRUE)

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./ITS_read_counts_at_each_step.csv", row.names = TRUE)


# IMPORT METADATA ####
meta <- read_csv("./em_community_metadata_clean.csv")
meta <- as.data.frame(meta)
row.names(meta) <- meta$SampleNumber
meta <- as_tibble(meta)

# reorder metadata to match seqtab
df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                 SampleNumber=row.names(seqtab.nochim))
df2 <- left_join(meta,df,by="SampleNumber")
df2 <- as.data.frame(df2)
row.names(df2) <- df2$SampleNumber
meta <- df2[row.names(seqtab.nochim),]
identical(row.names(meta),row.names(seqtab.nochim))


# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove singleton taxa - looks like it removed 14
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# remove newly singleton taxa - maybe removed 10 more? 
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./seqtab.nochim.clean.RDS")

seqtab.nochim <- readRDS("./seqtab.nochim.clean.RDS")

# Find and remove contaminants ####
# need to explain how to identify these in data set
meta$Control <- meta$Host_ID == "Blank"
meta <- meta[meta$SampleNumber %in% row.names(seqtab.nochim),] 

contams = isContaminant(seqtab.nochim, neg = meta$Control, normalize = TRUE)

table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./likely_contaminants.csv", row.names = TRUE)


# Inspect the control samples 
controls <- seqtab.nochim[194:198, ]

# Filter to only include columns that have value > 0
controls <- as.data.frame(controls)

controls <- controls[, colSums(controls > 0) > 0]

# transpose to make it easier to read 
controls <- t(controls)




# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Control == FALSE,]
meta = meta[meta$Control == FALSE,]


################### SAVE FILES FOR LATER USE ##########

#save modified meta file 
write.csv(meta, file = "./EM_meta_final.csv", row.names = TRUE)

#save seqtab with contaminants removed. Fully ready for taxonomy - 25,276 ASV's
saveRDS(seqtab.nochim,"./seqtab.nochim.final.RDS")

# ASSIGN TAXONOMY ####
#this step was performed in Kamiak due to resource constrains locally 

# downloaded the full database from UNITE, which is the choice for ITS sequences 
# this is stored in the 'data' file, so switch the working directory 
unite.ref <- "./UNITE_public_04.04.2024.fasta"
taxa <- assignTaxonomy(seqtab.nochim, "./UNITE_public_04.04.2024.fasta", multithread=4, tryRC=TRUE, verbose=TRUE)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./output/ITS_Taxonomy_from_dada2.RDS")

# add_species names - remove for fungal taxa, this looks through for exact matches for bacteria
#taxa <- addSpecies(taxa, "./taxonomy/rdp_species_assignment_16.fa.gz")

# Save completed taxonomy file
#saveRDS(taxa, file = "./output/ITS_Taxonomy_from_dada2_sp.RDS")

##############################################################

