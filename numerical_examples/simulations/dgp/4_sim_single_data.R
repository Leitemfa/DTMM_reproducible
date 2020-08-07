
library(MASS)

sim_single <- function(seed,
                          n1 = 40,
                          n2 = 30,
                          n3 = 20,
                          a = 3,
                          b = 1,
                          mu_1 = c(3, 1, a, b, 0),
                          mu_2 = c(2.434, 2.434, a, b, 0),
                          mu_3 = c(1, 3, a, b, 0),
                          sigma = 0.05,
                          sigma2 = 1,
                          save_data = FALSE,
                          file_dir){

  set.seed(seed)

  Sigma_1 <- matrix(c(sigma, 0, 0, 0, 0,
                     0, sigma, 0, 0, 0,
                     0, 0, sigma2, 0, 0,
                     0, 0, 0, sigma2, 0,
                     0, 0, 0, 0, sigma2), ncol = 5, nrow = 5, byrow = TRUE)

  Sigma_2 <- matrix(c(sigma, 0, 0, 0, 0,
                     0, sigma, 0, 0, 0,
                     0, 0, sigma2, 0, 0,
                     0, 0, 0, sigma2, 0,
                     0, 0, 0, 0, sigma2), ncol = 5, nrow = 5, byrow = TRUE)

  Sigma_3 <- matrix(c(sigma, 0, 0, 0, 0,
                     0, sigma, 0, 0, 0,
                     0, 0, sigma2, 0, 0,
                     0, 0, 0, sigma2, 0,
                     0, 0, 0, 0, sigma2), ncol = 5, nrow = 5, byrow = TRUE)

  x_1 <- mvrnorm(n1, mu_1, Sigma_1)
  x_2 <- mvrnorm(n2, mu_2, Sigma_2)
  x_3 <- mvrnorm(n3, mu_3, Sigma_3)

  N_1 <- rnbinom(n1, mu = 15000, size = 20)
  N_2 <- rnbinom(n2, mu = 15000, size = 20)
  N_3 <- rnbinom(n3, mu = 15000, size = 20)

  y_1 <- matrix(0, ncol = 6, nrow = n1)
  y_2 <- matrix(0, ncol = 6, nrow = n2)
  y_3 <- matrix(0, ncol = 6, nrow = n3)

  for(i in 1:n1){
    x <- x_1[i, ]
    a <- sum(exp(x)) + 1
    p <- c(exp(x), 1)/a
    y_1[i, ] <- rmultinom(1, size = N_1[i], prob = p)
  }


  for(i in 1:n2){
    x <- x_2[i, ]
    a <- sum(exp(x)) + 1
    p <- c(exp(x), 1)/a
    y_2[i, ] <- rmultinom(1, size = N_2[i], prob = p)
  }

  for(i in 1:n3){
    x <- x_3[i, ]
    a <- sum(exp(x)) + 1
    p <- c(exp(x), 1)/a
    y_3[i, ] <- rmultinom(1, size = N_3[i], prob = p)
  }

  y <- rbind(y_1, y_2, y_3)
  colnames(y) <- seq(1, 6)
  rownames(y) <- seq(1, n1 + n2 + n3)

  if(save_data == FALSE){
    return(y)
  }
  else{
    data <- y
    write.csv(data, file = file_dir)
    return(y)
  }
}
