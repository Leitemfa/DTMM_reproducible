
rm(list = ls())
setwd("~/data/AG")

library(data.table)
library(phyloseq)
library(ape)
library(dplyr)

AG_tree = read_tree_greengenes("97_otus.tree")        #-import the phylogenetric tree
otu = fread('ag_fecal_from_biom.txt', header = TRUE)  #-import the otu table


otu_matrix = do.call(rbind, otu)
otu_matrix = apply(otu_matrix,1, as.numeric)
row.names(otu_matrix) = otu_matrix[, "OTUID"]
otu_neat = otu_matrix[, -1]
rm(otu_matrix)


#################################################################################################
####----choose subset of samples

ag_fecal = fread("ag_fecal.txt")

#ag_fecal = ag_fecal[ag_fecal$DIABETES == "Diagnosed by a medical professional (doctor, physician assistant)"]
ag_fecal = ag_fecal[ag_fecal$IBD == "Diagnosed by a medical professional (doctor, physician assistant)"]

otu_neat = otu_neat[, intersect(colnames(otu_neat), ag_fecal$SampleID) ]

#################################################################################################
####----choose subset of otus based on some criterion

otu_rowsum = apply(otu_neat, 1, sum)
otu_colsum = apply(otu_neat, 2, sum)

num_otu_top = 75  ####----number of top otus to be used in the analysis
num_min_sample = 500

cutoff = sort(otu_rowsum, decreasing = TRUE)[num_otu_top]

otu_top = otu_neat[otu_rowsum >= cutoff, ]
otu_top = otu_top[, apply(otu_top, 2, sum) >= num_min_sample]

ag_fecal = ag_fecal[ag_fecal$SampleID %in% colnames(otu_top), ]

#################################################################################################
####---create a trimmed otu tree with only the otus selected above

diff = setdiff(AG_tree$tip.label, row.names(otu_top))
tree = drop.tip(AG_tree, tip=diff, trim.internal = TRUE, subtree = FALSE, root.edge = 0)

#################################################################################################
####----remove the ununsed variables

rm(otu, otu_neat, AG_tree, cutoff, diff, num_otu_top, otu_rowsum, otu_colsum)

#################################################################################################
####----save the data

save(otu_top, tree, ag_fecal, file="./AG_ibd.RData")
