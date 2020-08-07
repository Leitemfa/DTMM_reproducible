

rm(list = ls())
setwd("~/numerical_examples/validation")
load("~/data/Dethlefsen_Relman/abt.rda")

library(ape)
library(phyloseq)
library(DTMM)
library(DirichletMultinomial)


abt_5 <- tax_glom(abt, taxrank = rank_names(abt)[5])

####----F

otu <- otu_table(abt_5)[, 57:108]
tree <- phy_tree(abt_5)
plot(tree)
y_E <- t(otu)

res_time_E_5 <- DTMM(y_E, tree, tau_vec = 10 ^ seq(-1, 4, 0.5), theta_vec = seq(0.01, 0.99, 0.08),
                     init_gamma = rep(1,57), alpha = "default", mcmc_iter = 2500, select = TRUE)

save(res_time_E_5, file = "E_5.RData")
