# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Building and adding a phylogeny to the cleaned phyloseq object
# Author: Geoffrey Zahn
# Modified by: L. McKinley Nevins 
# April 2, 2025
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     vegan v 2.6.10
#                     phangorn v 2.12.1
#                     phyloseq v 1.48.0
#                     msa v 2.12.1
#                     ape v 5.8.1
#                     seqinr v 4.2.36
#                     ggtree v 3.12.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Perform multiple sequence alignment of all ASVs, build distance matrix,       # 
# construct and refine a phylogenetic tree, add the tree to the phyloseq object #
#           With larger data sets, this can be a long process...                #
# Further, proper phylogenetics is beyond the scope of this tutorial.           #
#################################################################################

# Packages and functions ####
BiocManager::install("msa")

BiocManager::install("ggtree")

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(ggtree); packageVersion("ggtree")

# Load phyloseq object produced from the subset of the funguild community 
# This includes just the AM taxa. Could use an earlier version of the phyloseq object prior 
# to filtering if we wanted to look at all of the taxa present in the community. 
ps <- readRDS("./FINAL/AM_phyloseq_final_2025.RDS")



ps2 <- 
ps %>% 
  subset_taxa(!is.na(tax_table(ps)[,1]))

# simplify OTU names
seqs <- rownames(tax_table(ps2))
names(seqs) <- paste0("OTU_",1:length(seqs)) # This propagates to the tip labels of the tree



# Multiple sequence alignment  ####
# This is with just ten iterations, could do more 
alignment <- msa::msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)

# save progress 
# saveRDS(alignment,"./output/18S_dna_alignment_muscle.RDS")

# Convert to phangorn format
phang.align = as.phyDat(alignment, type = "DNA")

## From here, this is trying something from this paper: 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC11117635/

# Comparing nucleotide substitution models 
modeltest <- modelTest(phang.align) 

# infering a tree using maximum likelihood method 
pml <- pml_bb(modeltest, method = "unrooted", start = NULL) 

# bootstrap analysis 
bstrees <- bootstrap.pml(pml, bs = 100, trees = TRUE, multicore = TRUE, mc.cores = 8) 

# started 

# plotting the ML tree with the bootstrap values 
plotBS(pml$tree, bstrees, type = "phylogram")


#####Get names of tree set to match with the phyloseq object 

# simplify OTU names of tree 
old_ids <- rownames(tax_table(ps2))  
new_ids <- paste0("OTU_", seq_along(old_ids))

# Create a named vector: names = old, values = new
id_map <- setNames(new_ids, old_ids)

#will need this later for identification purposes 
id_table <- as.data.frame(id_map)

# load the new ID's onto the tree 
pml$tree$tip.label <- new_ids

# load into the taxa names of the ps2 object 
taxa_names(ps2) <- id_map[taxa_names(ps2)]


## load tree into phyloseq object, now that names are matched 
ps_boot <- phyloseq(tax_table(tax_table(ps2)),
                    otu_table(otu_table(ps2)),
                    sample_data(sample_data(ps2)),
                    phy_tree(pml$tree))

saveRDS(ps_boot, "./clean_AM_w_tree_2025.RDS")

# This phyloseq object now has the maximum likelihood tree loaded in for plotting 

#################
# plotting with phyloseq 

# color by family 
family <- plot_tree(ps_boot, method = "sampledodge", color = "Family", label.tips = NULL, 
                    plot.margin = 0, ladderize = "left", base.spacing = 0.01)

family

family_treeonly <- plot_tree(ps_boot, method = "treeonly", color = "Family", 
                             plot.margin = 0, ladderize = "left", base.spacing = 0.01)

family_treeonly

# color by site 
site <- plot_tree(ps_boot, method = "sampledodge", color = "Site", label.tips = NULL, 
                  plot.margin = 0, ladderize = "left", base.spacing = 0.01)

site

####################
# try making a tree in RAxML - NG

# take previous alignment and save it in fasta format 

# not working, can revisit later 
write.dna(alignment, file = "AM_otus_aligned_2025.fasta", format = "fasta", nbcol = 264, colsep = "")





###################################################
## SAND BOX
###################################################
# old way, but more exploratory and not sufficient for a real phylogenetic tree 

# distance - maximum likelihood ####
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"./output/18S_ML_Distance.RDS")

# Initial neighbor-joining tree ####
treeNJ <- NJ(dm) # Note, tip order != sequence order

# save progress
saveRDS(treeNJ, "./output/18S_treeNJ.RDS")

# Estimate model parameters ####
fit = pml(treeNJ, data=phang.align)

#save
saveRDS(fit,"./output/18S_fit_treeNJ.RDS")


# add tree to phyloseq object ####
# running into some issues here with the names 
# looks like I edited the names further upstream, so now it's leading to 
# issues matching up down here 
ps2 <- phyloseq(tax_table(tax_table(ps2)),
                otu_table(otu_table(ps2)),
                sample_data(sample_data(ps2)),
                phy_tree(fit$tree))

plot_tree(ps2)



# simplify OTU names of tree 
old_ids <- rownames(tax_table(ps2))  
new_ids <- paste0("OTU_", seq_along(old_ids))

# Create a named vector: names = old, values = new
id_map <- setNames(new_ids, old_ids)

#will need this later for identification purposes 
id_table <- as.data.frame(id_map)

# load the new ID's onto the tree 
fit$tree$tip.label <- new_ids

plot(fit)


# Save updated phyloseq object with tree
#saveRDS(ps2, "./output/ps_not-cleaned_w_tree.RDS")

# plotting phyloseq tree 
ps2

AM_Tree <- phy_tree(ps2)

require(ggtree)

# ggtree running into a few issues because it doesn't have all of the same metadata 
# associated with it. For now would be better to just use Phyloseq plotting 
tree_plot <- ggtree(AM_Tree) + geom_point(aes(color=Family), size=3)

tree_plot


ggtree(AM_Tree, branch.length='none', layout='circular')

test <- ggtree(AM_Tree, layout='circular') + xlim(-0.05, NA)

test

#################
# plotting with phyloseq 

# color by family 
family <- plot_tree(ps2, method = "sampledodge", color = "Family", label.tips = NULL, 
          plot.margin = 0, ladderize = "left", base.spacing = 0.01)

family

# color by site 
site <- plot_tree(ps2, method = "sampledodge", color = "Site", label.tips = NULL, 
          plot.margin = 0, ladderize = "left", base.spacing = 0.005)

site

ps2@tax_table[,6]
