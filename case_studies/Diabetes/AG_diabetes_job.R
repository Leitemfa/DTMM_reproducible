
rm(list = ls())
library(ape)
library(DTMM)

load("~/data/AG/AG_diabetes.RData")
setwd("~/case_studies/Diabetes")

tree = reorder(tree, "postorder")
y = t(otu_top)

res_ag_diabetes = DTMM(y, tree, tau_vec = 10 ^ seq(-1, 4, 0.5), theta_vec = seq(0.01, 0.99, 0.08),
                       alpha = "default", mcmc_iter = 2500, select = TRUE)

save(res_ag_diabetes, file = "res_ag_diabetes.RData")

