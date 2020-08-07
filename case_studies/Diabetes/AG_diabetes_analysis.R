
rm(list = ls())

load("~/data/AG/AG_ibd.RData")
setwd("~/case_studies/Diabetes")
load("/res_ag_diabetes.RData")

library(DirichletMultinomial)
library(clusteval)
library(vegan)
library(ggplot2)
library(reshape2)
library(viridis)

B <- 2500 ##--burn in
Total <- 5000
n <- dim(res_ag_diabetes$post_c)[1]
Nsample <- dim(res_ag_diabetes$post_c)[2] - B


####----save posterior samples after burn-in
c <- res_ag_diabetes$post_c
gamma <- res_ag_diabetes$post_gamma
alpha <- res_ag_diabetes$post_alpha
lambda <- res_ag_diabetes$post_lambda


dmm_2 <- dmn(t(otu_top), k = 2)
dmm_cluster <- apply(dmm_2@group, 1, which.max)


####----Generate the treaceplots
par(mfrow = c(2, 2))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(1:Total, alpha, type = "l", ylab = expression(beta), xlab = "Iteration", main = "(a)")
abline(v = B, col = "blue", lty = 2)
abline(h = mean(alpha[(B+1):Total]), col = "red", lty = 2)
plot(1:Total, lambda, type = "l", ylab = expression(lambda), xlab = "Iteration", main = "(b)")
abline(v = B, col = "blue", lty = 2)
abline(h = mean(lambda[(B+1):Total]), col = "red", lty = 2)
plot(1:Total, apply(gamma, 2, sum), type = "l", ylab = expression(paste(Sigma,gamma(A))), xlab = "Iteration", main = "(c)")
abline(v = B, col = "blue", lty = 2)
abline(h = mean(apply(gamma, 2, sum)[(B+1):Total]), col = "red", lty = 2)

c1 <- rep(0, Total)
c2 <- rep(0, Total)
c3 <- rep(0, Total)
c4 <- rep(0, Total)
c5 <- rep(0, Total)
for(i in 1:Total){
  label_count <- sort(table(c[, i]), decreasing = TRUE)
  c1[i] <- label_count[1]/n
  c2[i] <- sum(label_count[1:2])/n
  c3[i] <- sum(label_count[1:3])/n
  c4[i] <- sum(label_count[1:4])/n
  c5[i] <- sum(label_count[1:5])/n
}

plot(1:Total, c5, type = "l", ylab = expression("ratio"), xlab = "Iteration", ylim = c(0.2, 1), main = "(d)")
lines(1:Total, c4, type = "l", xlab = "Iteration")
lines(1:Total, c3, type = "l", xlab = "Iteration")
lines(1:Total, c2, type = "l", xlab = "Iteration")
lines(1:Total, c1, type = "l", xlab = "Iteration")
abline(v = B, col = "blue", lty = 2)


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


partition_est[partition_est == 3] <- 2
partition_est[partition_est == 5] <- 3
partition_est = factor(partition_est)
l = length(levels(partition_est))
levels(partition_est) = seq(1, l)
re_order = order(partition_est)

y = t(otu_top)

y_bar = apply(y, 1, sum)
y_ratio = sweep(y, 1, y_bar, '/')
otu_order = rev(tree$tip.label)
y_sort = y_ratio[ , otu_order]
y_sort = y_sort[re_order, ]
rownames(y_sort) = factor(1:length(partition_est))

comelt <- melt((y_sort)^0.5)
comelt$cluster <- as.numeric(apply(as.matrix(comelt$Var1), 1, function(x){sort(partition_est)[x]}))
head(comelt)
comelt$Var1 = factor(comelt$Var1)
comelt$Var2 = factor(comelt$Var2, levels = rev(colnames(y_sort)))


####----Plot the heatmaps of samples grouped by their clustering labels.

p <- ggplot(data = comelt, mapping = aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_gradient(name = expression(tilde(y)[ij]^0.5),  low = "#FFFFFF", high = "black", na.value="white") +
  xlab(label = "Sample") + ylab(label = "OTU") +
  facet_grid(~ cluster, switch = "x", scales = "free_x", space = "free_x") +
  theme(strip.placement = "outside", axis.text.x =element_blank()) +theme_bw() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

p + geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 38.5, ymax = 39.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 54.5, ymax = 55.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 59.5, ymax = 60.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  #geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 25.5, ymax = 26.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 25.5, ymax = 26.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 31.5, ymax = 32.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 2), aes(xmin = 0.5, xmax = 42.5, ymin = 38.5, ymax = 39.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 2), aes(xmin = 0.5, xmax = 41.5, ymin = 54.5, ymax = 55.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 9.5, ymax = 10.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 18.5, ymax = 19.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 6.5, ymax = 7.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 2), aes(xmin = 0.5, xmax = 41.5, ymin = 18.5, ymax = 19.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE)



cocluster_mean_re = cocluster_mean[re_order, re_order]
comelt <- melt(cocluster_mean_re)
head(comelt)

####----Get the coclustering plot
ggplot(data = comelt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  theme_bw() +
  scale_fill_gradient2(low = "#FFFFFF", high = "black", name = expression(hat(Pi)[paste(i[1],i[2])])) +
  coord_fixed() +
  xlab(label = "Sample") + ylab(label = "Sample")+
  theme(strip.placement = "outside", axis.text.x =element_blank()) +theme_bw()



y = t(otu_top)

example_NMDS = metaMDS(y, k = 2)
r = partition_est

r = factor(r)
l = length(levels(r))
levels(r) = seq(1, l)

color = c("blue", "red", "darkgoldenrod1")
pch = c(15,16,17)

mycol = rep(0, length(r))
mypch = rep(0, length(r))
for(i in 1:length(mycol)){
  mycol[i] = color[as.numeric(r[i])]
  mypch[i] = pch[as.numeric(r[i])]
}

####----Generate the 2D NMDS plot for DTMM
plot(example_NMDS, display = "sites", type = "n", main = "DTMM" , xlim = c(-0.8, 1.1))
points(example_NMDS, disp="sites", pch = mypch, col=mycol, bg = mycol, cex=0.8)
legend("topright", legend = c("cluster 1", "cluster 2", "cluster 3"), col = color,
       pch = pch, cex = 0.8)


dmm_2 <- dmn(t(otu_top), k = 2)
dmm_cluster <- apply(dmm_2@group, 1, which.max)
r = dmm_cluster
r = factor(r)
l = length(levels(r))
levels(r) = seq(1, l)

color = c("blue", "red", "darkgoldenrod1")
pch = c(15,16)
mycol = rep(0, length(r))
mypch = rep(0, length(r))
for(i in 1:length(mycol)){
  mycol[i] = color[as.numeric(r[i])]
  mypch[i] = pch[as.numeric(r[i])]
}

####----Generate the 2D NMDS plot for DMM
plot(example_NMDS, display = "sites", type = "n", main = "DMM", xlim = c(-0.8, 1.1))
points(example_NMDS, disp="sites", pch = mypch, col = mycol, bg = mycol, cex=0.8)
legend("topright", legend = c("cluster A", "cluster B"), col = c("blue", "red"),
       pch = c(15, 16), cex = 0.8)


####
fit <- mclapply(1:7, dmn, count=t(otu_top), verbose=TRUE)
lplc <- sapply(fit, laplace)

####----Generate the fitting plot for DMM
plot(1:7, lplc, type="l", xlab="Number of Dirichlet Components", ylab= "Model Fit", main = "Diabetes")
points(1:7, lplc, pch = 19)


gamma_mean = apply(gamma[,(B+1):Total], 1, mean)

####----Plot the posterior node selection probability for each node
par(mai=c(0.5, 0.4, 0.6 , 0))

col_Pal = viridis_pal(option = "A", direction = -1)
node_col = col_Pal(500)[as.numeric(cut(c(gamma_mean, 0, 1), breaks = 500)) ]
plot(tree, show.tip.label = TRUE, use.edge.length = FALSE, show.node.label=FALSE,
     cex = 0.6)
nodelabels(text=format(round(gamma_mean, digits=2), nsmall=2), cex=0.2, col = node_col, bg=node_col, frame="circle")




####---Compute the importance for each OTU in the clustering procedure

tax = read.table("~/data/AG/97_otu_taxonomy.txt", sep = ";", header = FALSE, fill = TRUE, colClasses = "character")

otu_names = rownames(otu_top)

text = tax[, 1]
elems = unlist( strsplit(tax[, 1], "\t" ) )
first = matrix( elems , ncol = 2 , byrow = TRUE )
tax = tax[, -1]
tax = cbind(first, tax)
rownames(tax) = tax[, 1]
tax = tax[intersect(otu_names, tax[ , 1]),  ]
tax = tax[ , -1]
colnames(tax) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax = as.matrix(tax)
rm(first, elems, text)

####

y = t(otu_top)
ybar = apply(y, 1, sum)
y = sweep(y, 1, ybar, '/')

ratio = rep(0, 75)

partition_est_s <- as.numeric(partition_est)
partition_est_s[partition_est_s != 1] <- 5
partition_est_s <- factor(partition_est_s)

for(i in 1:75){
  otu_split = split(y[, i], partition_est_s)
  within_var = sum(unlist(lapply(otu_split, var)) * (table(partition_est_s) - 1), na.rm = TRUE)
  between_var = sum( (unlist(lapply(otu_split, mean)) - mean(y[, i]))^2 * table(partition_est_s))
  total_var = num = var(y[, i]) * (sum(table(partition_est_s)) - 1)
  print(c(colnames(y)[i], within_var, between_var, total_var))
  ratio[i] = between_var/within_var
}

data = data.frame(rownames(otu_top), ratio)
data = data[order(data$ratio, decreasing = TRUE)[1:10], ]
data


table(partition_est)

p + geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 38.5, ymax = 39.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
    geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 54.5, ymax = 55.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
    geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 59.5, ymax = 60.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
    #geom_rect(data = data.frame(cluster = 1), aes(xmin = 0.5, xmax = 41.5, ymin = 25.5, ymax = 26.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
    geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 25.5, ymax = 26.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 31.5, ymax = 32.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
    geom_rect(data = data.frame(cluster = 2), aes(xmin = 0.5, xmax = 42.5, ymin = 38.5, ymax = 39.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
    geom_rect(data = data.frame(cluster = 2), aes(xmin = 0.5, xmax = 41.5, ymin = 54.5, ymax = 55.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 9.5, ymax = 10.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 18.5, ymax = 19.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 3), aes(xmin = 0.5, xmax = 23.5, ymin = 6.5, ymax = 7.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE) +
  geom_rect(data = data.frame(cluster = 2), aes(xmin = 0.5, xmax = 41.5, ymin = 18.5, ymax = 19.5), alpha = 0, color = "black", size = 0.7, inherit.aes = FALSE)



sub_tax = tax[as.character(data[, 1]),]


show_table = round(t(C[,as.character(data[, 1])] * 100), digits = 2)
show_table = cbind(round(data[,2],digits = 2), show_table)

tax[row.names(show_table), 5]

table(tax[,5])

row.names(tax[tax[,6] == " g__Ruminococcus", ])
row.names(tax[tax[,6] == " g__Bacteroides", ])
row.names(tax[tax[,6] == " g__Prevotella", ])



tax = tax[otu_order,]



centroid = HDTM_centroid(t(otu_top), tree, c = partition_est, gamma = gamma[,B + which.min(least_square) ])


c1 = centroid$centroid[[1]]
c2 = centroid$centroid[[2]]
c3 = centroid$centroid[[3]]



C = t(cbind(c1,c2,c3))
colnames(C) = colnames(y)
C = C[, otu_order]
rownames(C) = c(1,2,3)

comelt <- melt((C)^0.5)
comelt$cluster <- as.numeric(rep(c(1,2,3), 75))
head(comelt)
comelt$Var1 = factor(comelt$Var1)
comelt$Var2 = factor(comelt$Var2, levels = rev(colnames(y_sort)))

####----Plot the estimated clustering centroids
ggplot(data = comelt, mapping = aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_viridis(option="magma", direction = -1, na.value="white", limits = c(0,1), name = expression(bar(p)[k]^0.5)) +
  xlab(label = "Cluster") + ylab(label = NULL) +
  facet_grid(~ cluster, switch = "x", scales = "free_x", space = "free_x") +
  theme(strip.placement = "outside", axis.text.x =element_blank()) +theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank())


y = t(otu_top)
ybar = apply(y, 1, sum)
y = sweep(y, 1, ybar, '/')

y_family <- y^(1)
rownames(y_family) <- seq(1, dim(y)[1])
tax_mod <- tax
tax_mod[tax_mod[,5] == " f__Rikenellaceae", 6] <- " f__Rikenellaceae"
colnames(y_family) <- tax_mod[colnames(y), 6]

y_family <- t(rowsum(t(y_family), group = colnames(y_family), na.rm = T))

keeps <- c( " g__Bacteroides", " g__[Ruminococcus]", " g__Prevotella", " g__Faecalibacterium", " g__Blautia", " f__Rikenellaceae")

y_family <- y_family[, keeps]
y_family <- cbind(y_family, y[,"173876"])
colnames(y_family)[7] <- "173876"
y_family <- y_family[order(y_family[,1]),]

comelt <- melt(y_family)
comelt$cluster <- as.numeric(apply(as.matrix(comelt$Var1), 1, function(x){partition_est[x]}))
head(comelt)
comelt$Var1 = factor(comelt$Var1, levels = row.names(y_family))
comelt$Var2 = factor(comelt$Var2, levels = names(sort(apply(y_family,2,sum))))

levels(comelt$Var2)

label_genus <- c("Ruminococcus", "173876", "Rikenellaceae*", "Blautia", "Prevotella", "Faecalibacterium", "Bacteroides")

brewer.pal(8, "Dark2")

####----Plot the relative abundance for representative Genus.
q <- ggplot(data = comelt, mapping = aes(x=Var1, y=value, fill=Var2)) + geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#D95F02","#7570B3", "#66A61E", "#E7298A", "#E6AB02", "#A6761D", "#666666"), name = "Genus", label = label_genus ) +
  xlab(label = "Cluster") + ylab(label = "Ratio") +
  facet_grid(~ cluster, switch = "x", scales = "free_x", space = "free_x") +
  theme(strip.placement = "outside", axis.text.x =element_blank()) +theme_bw() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

alpha_div <- data.frame(Shannon = diversity(y, index = "shannon"), cluster = partition_est)

####----Plot the Shannon index for each cluster
ggplot(alpha_div, aes( cluster, Shannon)) + geom_boxplot()+theme_bw()
















