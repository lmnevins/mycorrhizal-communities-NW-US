# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Processing raw amplicon reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     dada2 v 1.18.0
#                     decontam v 1.6.0
#                     phyloseq v 1.30.0
#                     purrr v 1.0.2
#                     Biostrings v 2.66.0
#                     patchwork v 1.0.1
# -----------------------------------------------------------------------------#

##################################################################
#### This script calls cutadapt to remove flanking ITS regions ####
#### You must have cutadapt installed on your system           ####
#### and present in your PATH. See cutadapt installation       ####
#### documents for instructions.                               ####
##################################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")


Sys.setenv(PATH="/usr/local/Caskroom/miniforge/base/envs/Cutadapt/bin:/opt/anaconda3/condabin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/Library/Apple/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin")

# PARSE FILE PATHS ####

# File parsing
path <- "/Users/mckinleynevins/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Data/240223_Nevins__ITS_demultiplxed/" # CHANGE to the subdirectory containing your demultiplexed fastq files
fnFs <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_001.fastq.gz")) # make pattern match your FWD reads
fnRs <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R2_001.fastq.gz")) # make pattern match your REV reads

  # PRIMER SITES WERE ALREADY REMOVED BY SEQUENCING FACILITY, BUT DOUBLECHECKING THAT THEY AREN'T PRESENT ####

######################################################################################################

# Here, you should supply the primer sequences used during PCR
FWD <- "CTTGGTCATTTAGAGGAAGTA" # Sequence of FWD primer
REV <- "GCTGCGTTCTTCATCGATGC"  # Sequence of REV primer
######################################################################################################

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients


# Prefilter to remove reads with ambiguous (N) bases ####
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterEd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
#TIME WARNING ON THIS STEP
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# We are now ready to count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations. Identifying and counting the primers on one 
# set of paired end FASTQ files is sufficient, assuming all the files were created using the same 
# library preparation, so we’ll just process the first sample.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#there are some hits on the reverse compliment - 4,911 for FWD.ReverseReads and 8,017 for REV.ForwardReads

# Run cutadapt ####
# If the following command returns an error, you do not have cutadapt installed correctly

system2("cutadapt", args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#TIME WARN - this takes a good bit of time to run and spits out a ton of output into the console 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# sanity check
# This should show no occurences in any of the orientations now
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#all gone now
#these end up in a cutadapt file now, which should be used for downstream applications
