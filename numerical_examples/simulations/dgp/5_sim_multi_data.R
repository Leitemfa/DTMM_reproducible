
sim_multi <- function(seed,
                      n1 = 40,
                      n2 = 30,
                      n3 = 20,
                      a = 3,
                      b = 3,
                      sigma = 0.05,
                      mu_1 = c(a, b, 3.5, 2.5, 3),
                      mu_2 = c(a, b, 2.5, 3, 3.5),
                      mu_3 = c(a, b, 3, 3.5, 2.5),
                      save_data = FALSE,
                      file_dir){

  set.seed(seed)

  Sigma_1 <- matrix(c(1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, sigma, 0, 0,
                     0, 0, 0, sigma, 0,
                     0, 0, 0, 0, sigma), ncol = 5, nrow = 5, byrow = TRUE)

  Sigma_2 <- matrix(c(1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, sigma, 0, 0,
                     0, 0, 0, sigma, 0,
                     0, 0, 0, 0, sigma), ncol = 5, nrow = 5, byrow = TRUE)

  Sigma_3 <- matrix(c(1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, sigma, 0, 0,
                     0, 0, 0, sigma, 0,
                     0, 0, 0, 0, sigma), ncol = 5, nrow = 5, byrow = TRUE)

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
