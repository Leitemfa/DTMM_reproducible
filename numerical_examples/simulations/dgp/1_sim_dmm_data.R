
library(gtools)

sim_dmm = function(seed,
                    n1 = 40,
                    n2 = 30,
                    n3 = 20,
                    alpha1 = c(2,2,3,2,2,2),
                    alpha2 = c(2,5,1,2,1,2),
                    alpha3 = c(4,3,2,1,2,1),
                    save_data = FALSE,
                    file_dir){

  set.seed(seed)

  N_1 = rnbinom(n1, mu = 15000, size = 20)
  N_2 = rnbinom(n2, mu = 15000, size = 20)
  N_3 = rnbinom(n3, mu = 15000, size = 20)

  y_1 = matrix(0, ncol = 6, nrow = n1)
  y_2 = matrix(0, ncol = 6, nrow = n2)
  y_3 = matrix(0, ncol = 6, nrow = n3)

  for(i in 1:n1){
    p <- rdirichlet(1, alpha1)
    y_1[i, ] = rmultinom(1, N_1[i], p)
  }


  for(i in 1:n2){
    p <- rdirichlet(1, alpha2)
    y_2[i, ] = rmultinom(1, N_2[i], p)
  }

  for(i in 1:n3){
    p <- rdirichlet(1, alpha3)
    y_3[i, ] = rmultinom(1, N_3[i], p)
  }

  y = rbind(y_1, y_2, y_3)
  colnames(y) <- seq(1, 6)
  rownames(y) <- seq(1, n1 + n2 + n3)

  if(save_data == FALSE){
    return(y)
  }
  else{
    data = y
    write.csv(data, file = file_dir)
    return(y)
  }
}

