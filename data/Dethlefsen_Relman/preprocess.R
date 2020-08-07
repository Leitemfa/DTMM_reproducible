

rm(list = ls())

library(ape)
library(phyloseq)
library(DirichletMultinomial)

#setwd("/data/Dethlefsen_Relman") ##--set the working directory to the data folder first
load("abt.rda")
abt_5 <- tax_glom(abt, taxrank = rank_names(abt)[5])

tree <- phy_tree(abt_5)

####----D

otu <- otu_table(abt_5)[, 1:56]
y_D <- t(otu)

####----E
otu <- otu_table(abt_5)[, 57:108]
y_E <- t(otu)

####----F
otu <- otu_table(abt_5)[,109:162]
y_F <- t(otu)

