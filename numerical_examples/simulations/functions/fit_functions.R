
library("DTTM")
library(ape)
library(cluster)
library(vegan)
library(kernlab)
library(DirichletMultinomial)

#install.packages("ape")
#install.packages("cluster")
#install.packages("vegan")
#install.packages("kernlab")
#install.packages("DirichletMultinomial")

####----fit DTMM
fit_hdtm <- function(y, tree, iter = 2){

  res_hdtm <- DTMM(y, tree, tau_vec = 10 ^ seq(-1, 4, 0.5), theta_vec = seq(0.01, 0.99, 0.08),
                   alpha = "default", mcmc_iter = iter, select = TRUE)

  B <- 1000
  Total <- 2000
  c <- res_hdtm$post_c[ , (B+1):Total]
  Nsample <- Total - B

  get_cocluster <- function(post_sample){
    ns <- dim(post_sample)[2]
    n <- dim(post_sample)[1]
    cocluster <- array(0, c(ns, n, n))
    for(i in 1:ns){
      for(j in 1:n){
        for(k in 1:n){
          cocluster[i, j, k] <- post_sample[j, i] == post_sample[k, i]
        }
      }
    }
    return(cocluster)
  }

  cocluster <- get_cocluster(c)
  cocluster_mean <- apply(cocluster, c(2, 3), mean)

  for(i in 1:dim(cocluster_mean)[1]){
    cocluster_mean[i, i] <- NA
  }

  least_square <- rep(0, Nsample)

  for(i in 1:Nsample){
    least_square[i] <- sum((cocluster[i, , ] - cocluster_mean) ^ 2, na.rm = TRUE)
  }

  partition_est <- c[, which.min(least_square)]

  cluster_hdtm <- partition_est
  return(cluster_hdtm)
}


####----fit DMM
fit_dmm <- function(y){

  dmm_fit <- mclapply(1:7, function(x){dmn(y, k = x)})
  num_of_cluster <- which.min(sapply(dmm_fit, laplace))
  cluster_dmm <- apply(dmm_fit[[num_of_cluster]]@group, 1, which.max)
  return(cluster_dmm)
}


####----fit Kmeans
fit_kmeans <- function(y, K = 3){
  cluster_kmeans <- kmeans(y, K)$cluster
  return(cluster_kmeans)
}


####----fit PAM with BC distance
fit_pam <- function(y, K = 3){
  y_rel <- sweep(y, 1, apply(y, 1, sum), "/")
  dy <- vegdist(y_rel)
  cluster_pam <- pam(dy, K)$clustering
  return(cluster_pam)
}


####----fit hierarchical clustering with BC distance
fit_hclust <- function(y, K = 3){
  y_rel <- sweep(y, 1, apply(y, 1, sum), "/")
  dy <- vegdist(y_rel)
  c <- hclust(dy)
  cluster_hclust <- cutree(c, K)
  return(cluster_hclust)
}


####----fit Dspectral clustering with chosen distance
fit_spec <- function(y, K = 3, distance = "default"){

  bc_diss <- function(x1, x2){
    return(vegdist(rbind(x1, x2)))
  }
  y_rel <- sweep(y, 1, apply(y, 1, sum), "/")

  if(distance != "bc"){
    c <- specc(y_rel, centers = K)
    cluster_spec <- c@.Data
  }else{
    c <- specc(y_rel, centers = K, kernal = "bc_diss")
    cluster_spec <- c@.Data
  }
  return(cluster_spec)
}


