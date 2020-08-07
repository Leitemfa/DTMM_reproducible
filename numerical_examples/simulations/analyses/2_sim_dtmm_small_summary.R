
rm(list = ls())

setwd("./simulations")
source("/functions/compute_similarity.R")

result_n <- array(0, c(100, 7, 90))
result_s <- array(0, c(100, 7, 90))
result_m <- array(0, c(100, 7, 90))
result_l <- array(0, c(100, 7, 90))

for(i in 0:99){
  file <- paste("results/res_dtmm_small_", i, ".RData", sep = "")
  load(file)
  result_n[i + 1, , ] <- res_dtmm$clusters_n
  result_s[i + 1, , ] <- res_dtmm$clusters_s
  result_m[i + 1, , ] <- res_dtmm$clusters_m
  result_l[i + 1, , ] <- res_dtmm$clusters_l
}

jaccard_n <- compute_similarity(result = result_n)
jaccard_s <- compute_similarity(result = result_s)
jaccard_m <- compute_similarity(result = result_m)
jaccard_l <- compute_similarity(result = result_l)

apply(jaccard_n, 1, mean)
apply(jaccard_s, 1, mean)
apply(jaccard_m, 1, mean)
apply(jaccard_l, 1, mean)

round(sqrt(apply((jaccard_n - 1)^2, 1, mean)), digits = 2)
round(sqrt(apply((jaccard_s - 1)^2, 1, mean)), digits = 2)
round(sqrt(apply((jaccard_m - 1)^2, 1, mean)), digits = 2)
round(sqrt(apply((jaccard_l - 1)^2, 1, mean)), digits = 2)

par(mfrow = c(2,2))
boxplot(t(jaccard_n), ylim = c(0, 1))
boxplot(t(jaccard_s), ylim = c(0, 1))
boxplot(t(jaccard_m), ylim = c(0, 1))
boxplot(t(jaccard_l), ylim = c(0, 1))


library(ggplot2)

jaccard_index1 <- t(jaccard_s)[ ,1:6]
jaccard_index2 <- t(jaccard_m)[ ,1:6]
jaccard_index3 <- t(jaccard_l)[ ,1:6]

jaccard_index = rbind(jaccard_index1, jaccard_index2, jaccard_index3)
method <- factor(c(rep("DTMM", 300), rep("DMM", 300), rep("K-means", 300), rep("PAM", 300), rep("Hclust", 300), rep("Spec", 300)),
                 levels = c("DTMM", "DMM", "K-means", "PAM", "Hclust", "Spec"))
signal <- factor(rep(c(rep("W", 100), rep("M", 100), rep("S", 100)), 6),
                 levels = c("W", "M", "S"))
jac <- c(jaccard_index[,1], jaccard_index[,2], jaccard_index[,3], jaccard_index[,4], jaccard_index[,5], jaccard_index[,6])
data <- data.frame(Jaccard = jac, Method = method, Signal = signal)

p <- ggplot(data, aes(x = Signal, y = Jaccard, fill = Method)) + facet_grid(~ Signal, switch = "x", scales = "free_x", space = "free_x") +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette="RdBu") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits=c(0.15,1)) + ggtitle("I. DT")

#ggsave("sec_3_60_jac_dtmm.eps", plot = p, width = 6, height = 5.4, path = "~/plots")



