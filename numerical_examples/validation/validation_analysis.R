
rm(list = ls())
library(ggplot2)
library(reshape2)

load("~/data/Dethlefsen_Relman/abt.rda")
setwd("~/numerical_examples/validation")
library(ape)
library(phyloseq)
library(DirichletMultinomial)

load("~/D_5.RData")
load("~/E_5.RData")
load("~/F_5.RData")

tax_table(abt)[tax_table(abt)[,4] == "Bacteroidaceae_Bacteroides", 5] <- "Bacteroides"
abt_5 <- tax_glom(abt, taxrank = rank_names(abt)[5])
tree <- phy_tree(abt_5)
plot(tree)

####----D

otu <- otu_table(abt_5)[,1:56]
y_D <- t(otu)

####----E

otu <- otu_table(abt_5)[,57:108]
y_E <- t(otu)


####----F

otu <- otu_table(abt_5)[,109:162]
y_F <- t(otu)

###################################################################################################


y <- y_D

fit <- mclapply(1:7, dmn, count=y, verbose=TRUE)
lplc <- sapply(fit, laplace)

plot(1:7, lplc, type="l", xlab="Number of Dirichlet Components", ylab= "Model Fit", main = "IBD")
points(1:7, lplc, pch = 19)

dmm_temp <- dmn(y, k = 2)
dmm_D <- apply(dmm_temp@group, 1, which.max)
dmm_D

####

y <- y_E

fit <- mclapply(1:7, dmn, count=y, verbose=TRUE)
lplc <- sapply(fit, laplace)

plot(1:7, lplc, type="l", xlab="Number of Dirichlet Components", ylab= "Model Fit", main = "IBD")
points(1:7, lplc, pch = 19)

dmm_temp <- dmn(y, k = 2)
dmm_E <- apply(dmm_temp@group, 1, which.max)
dmm_E

####


y <- y_F

fit <- mclapply(1:7, dmn, count=y, verbose=TRUE)
lplc <- sapply(fit, laplace)

plot(1:7, lplc, type="l", xlab="Number of Dirichlet Components", ylab= "Model Fit", main = "IBD")
points(1:7, lplc, pch = 19)

dmm_temp <- dmn(y, k = 2)
dmm_F <- apply(dmm_temp@group, 1, which.max)
dmm_F



B <- 1250
Total <- 2500
Nsample <- Total - B

####---get co-cluster matrix
get_cocluster <- function(post_sample){
  ns = dim(post_sample)[2]
  n = dim(post_sample)[1]
  cocluster = array(0, c(ns, n, n))
  for(i in 1:ns){
    for(j in 1:n){
      for(k in 1:n){
        cocluster[i, j, k] = post_sample[j, i] == post_sample[k, i]
      }
    }
  }
  return(cocluster)
}

####
c <- res_time_F_5$post_c

cocluster = get_cocluster(c[, (B+1):Total])
cocluster_mean = apply(cocluster, c(2, 3), mean)

for(i in 1:dim(cocluster_mean)[1]){
  cocluster_mean[i,i] = NA
}

least_square = rep(0, Nsample)

for(i in 1:Nsample){
  least_square[i] = sum((cocluster[i, , ] - cocluster_mean) ^ 2, na.rm = TRUE)
}

partition_est = c[, B + which.min(least_square)]

y <- y_F
y_bar = apply(y, 1, sum)
y_ratio = sweep(y, 1, y_bar, '/')
y_sort <- y_ratio

comelt <- melt((y_sort)^0.5)
comelt$cluster <- as.numeric(apply(as.matrix(comelt$Var1), 1, function(x){sort(dmm_cluster)[x]}))
head(comelt)
comelt$Var1 = factor(comelt$Var1)
comelt$Var2 = factor(comelt$Var2, levels = rev(colnames(y_sort)))

col <- c("firebrick", "steelblue", "darkgoldenrod2", "darkorchid", "yellow")
a <- partition_est
for(i in 1:length(partition_est)){
  a[i] <- col[partition_est[i]]
  #a[i] <- col[dmm_F[i]]
}


p_f1 <- ggplot(data = comelt, mapping = aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_gradient(name = expression(tilde(y)[ij]^0.5),  low = "#FFFFFF", high = "black", na.value="white") +
  xlab(label = "Sample") + ylab(label = "OTU") + theme_bw() + theme(axis.text.y=element_text(size=5), axis.text.x=element_text(size=8)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))

pdtmm_f <- p_f1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a)) +
  geom_vline(xintercept = c(11.5, 16.5, 40.5, 44.5), col = "darkblue") +
  geom_vline(xintercept = c(23.5, 51.5), col = "darkred") +  ggtitle("Subject F: DTMM")  +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 0, xmax = 11.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 11.5, xmax = 16.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 16.5, xmax = 23.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 23.5, xmax = 40.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 40.5, xmax = 44.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 44.5, xmax = 51.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 51.5, xmax = 55, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  annotate(geom = "text", x = 1, y = 61.5, label = "Pre", hjust = 0) +
  annotate(geom = "text", x = 12, y = 61.5, label = "CP1", hjust = 0) +
  annotate(geom = "text", x = 17, y = 61.5, label = "WPC1", hjust = 0) +
  annotate(geom = "text", x = 24, y = 61.5, label = "Interim", hjust = 0) +
  annotate(geom = "text", x = 41, y = 61.5, label = "CP2", hjust = 0) +
  annotate(geom = "text", x = 45, y = 61.5, label = "WPC2", hjust = 0) +
  annotate(geom = "text", x = 52, y = 61.5, label = "Post", hjust = 0) +
  theme(plot.title = element_text(hjust = 0.5))



col <- c("firebrick", "steelblue", "darkgoldenrod2", "black", "yellow")
a <- partition_est
for(i in 1:length(partition_est)){
  #a[i] <- col[partition_est[i]]
  a[i] <- col[dmm_F[i]]
}

pdmm_f <- p_f1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a)) +
  geom_vline(xintercept = c(11.5, 16.5, 40.5, 44.5), col = "darkblue") +
  geom_vline(xintercept = c(23.5, 51.5), col = "darkred") +  ggtitle("Subject F: DMM")  +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 0, xmax = 11.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 11.5, xmax = 16.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 16.5, xmax = 23.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 23.5, xmax = 40.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 40.5, xmax = 44.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 44.5, xmax = 51.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 51.5, xmax = 55, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  annotate(geom = "text", x = 1, y = 61.5, label = "Pre", hjust = 0) +
  annotate(geom = "text", x = 12, y = 61.5, label = "CP1", hjust = 0) +
  annotate(geom = "text", x = 17, y = 61.5, label = "WPC1", hjust = 0) +
  annotate(geom = "text", x = 24, y = 61.5, label = "Interim", hjust = 0) +
  annotate(geom = "text", x = 41, y = 61.5, label = "CP2", hjust = 0) +
  annotate(geom = "text", x = 45, y = 61.5, label = "WPC2", hjust = 0) +
  annotate(geom = "text", x = 52, y = 61.5, label = "Post", hjust = 0) +
  theme(plot.title = element_text(hjust = 0.5))

pdmm_f



####


c <- res_time_D_5$post_c

cocluster = get_cocluster(c[, (B+1):Total])
cocluster_mean = apply(cocluster, c(2, 3), mean)

for(i in 1:dim(cocluster_mean)[1]){
  cocluster_mean[i,i] = NA
}

least_square = rep(0, Nsample)

for(i in 1:Nsample){
  least_square[i] = sum((cocluster[i, , ] - cocluster_mean) ^ 2, na.rm = TRUE)
}

partition_est <- c[, B + which.min(least_square)]
partition_est <- factor(partition_est)
levels(partition_est) <- c(1,2,3)
partition_est <- as.numeric(partition_est)

y <- y_D
y_bar = apply(y, 1, sum)
y_ratio = sweep(y, 1, y_bar, '/')
y_sort <- y_ratio

get_label_D <- function(x){
  if(x <= 11){return(1)}
  else if(x <= 16){return(2)}
  else if(x <= 23){return(3)}
  else if(x <= 39){return(4)}
  else if(x <= 44){return(5)}
  else if(x <= 51){return(6)}
  else{return(7)}
}


comelt <- melt((y_sort)^0.5)
comelt$cluster <- as.numeric(apply(as.matrix(as.numeric(comelt$Var1)), 1, get_label_D))
head(comelt)
comelt$Var1 = factor(comelt$Var1)
comelt$Var2 = factor(comelt$Var2, levels = rev(colnames(y_sort)))

col <- c("darkgoldenrod2", "firebrick", "steelblue", "black", "yellow")
a <- partition_est
for(i in 1:length(partition_est)){
  a[i] <- col[partition_est[i]]
  #a[i] <- col[dmm_F[i]]
}

p_d1 <- ggplot(data = comelt, mapping = aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_gradient(name = expression(tilde(y)[ij]^0.5),  low = "#FFFFFF", high = "black", na.value="white") +
  xlab(label = "Sample") + ylab(label = "OTU") + theme_bw() + theme(axis.text.y=element_text(size=5), axis.text.x=element_text(size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))


pdtmm_d <- p_d1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a)) +
  geom_vline(xintercept = c(11.5, 16.5, 39.5, 44.5), col = "darkblue") +
  geom_vline(xintercept = c(23.5, 51.5), col = "darkred") +  ggtitle("Subject D: DTMM")  +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 0, xmax = 11.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 11.5, xmax = 16.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 16.5, xmax = 23.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 23.5, xmax = 39.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 39.5, xmax = 44.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 44.5, xmax = 51.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 51.5, xmax = 57, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  annotate(geom = "text", x = 1, y = 61.5, label = "Pre", hjust = 0) +
  annotate(geom = "text", x = 12, y = 61.5, label = "CP1", hjust = 0) +
  annotate(geom = "text", x = 17, y = 61.5, label = "WPC1", hjust = 0) +
  annotate(geom = "text", x = 24, y = 61.5, label = "Interim", hjust = 0) +
  annotate(geom = "text", x = 40, y = 61.5, label = "CP2", hjust = 0) +
  annotate(geom = "text", x = 45, y = 61.5, label = "WPC2", hjust = 0) +
  annotate(geom = "text", x = 52, y = 61.5, label = "Post", hjust = 0) +
  theme(plot.title = element_text(hjust = 0.5))



col <- c("steelblue", "firebrick", "darkgoldenrod2", "black", "yellow")
a <- partition_est
for(i in 1:length(partition_est)){
  #a[i] <- col[partition_est[i]]
  a[i] <- col[dmm_D[i]]
}

#pdmm_d <- p_d1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a)) +
#  geom_vline(xintercept = c(11.5, 16.5, 39.5, 44.5), col = "darkblue") +
#  geom_vline(xintercept = c(23.5, 51.5), col = "darkred") +  ggtitle("Subject D: DMM")  + theme(plot.title = element_text(hjust = 0.5))

pdmm_d <- p_d1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a)) +
  geom_vline(xintercept = c(11.5, 16.5, 39.5, 44.5), col = "darkblue") +
  geom_vline(xintercept = c(23.5, 51.5), col = "darkred") +  ggtitle("Subject D: DMM")  +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 0, xmax = 11.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 11.5, xmax = 16.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 16.5, xmax = 23.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 23.5, xmax = 39.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 39.5, xmax = 44.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 44.5, xmax = 51.5, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 51.5, xmax = 57, ymin = 59.5, ymax = 63.5), alpha = 1, fill = "azure2", color = "black",  inherit.aes = FALSE) +
  annotate(geom = "text", x = 1, y = 61.5, label = "Pre", hjust = 0) +
  annotate(geom = "text", x = 12, y = 61.5, label = "CP1", hjust = 0) +
  annotate(geom = "text", x = 17, y = 61.5, label = "WPC1", hjust = 0) +
  annotate(geom = "text", x = 24, y = 61.5, label = "Interim", hjust = 0) +
  annotate(geom = "text", x = 40, y = 61.5, label = "CP2", hjust = 0) +
  annotate(geom = "text", x = 45, y = 61.5, label = "WPC2", hjust = 0) +
  annotate(geom = "text", x = 52, y = 61.5, label = "Post", hjust = 0) +
  theme(plot.title = element_text(hjust = 0.5))



