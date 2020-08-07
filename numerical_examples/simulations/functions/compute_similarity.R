
library(clusteval)

####--result is a T * 7 * n array that saves all the simulation clustering results

compute_similarity <- function(result, measure = "jaccard"){

  num_sim <- dim(result)[1]
  num_method <- dim(result)[2]
  similarity <- matrix(0, ncol = num_sim, nrow = num_method)

  for(t in 1:num_sim){
    true_index <- result[t, num_method, ]
    similarity[ , t] <- apply(result[t, , ], 1, function(x){cluster_similarity(x, true_index, similarity = measure)})
  }
  return(similarity)
}

