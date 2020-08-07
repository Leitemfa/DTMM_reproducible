

sim_dt <- function(seed,
                   n,
                   beta1, beta2, beta3, beta4,
                   alpha1, alpha2, alpha3, alpha4, alpha5, alpha6){

  set.seed(seed)
  p <- matrix(0, ncol = 6, nrow = n)

  for(i in 1:n){
    b_left = rbeta(1, beta1, beta2)
    b_left_left = rbeta(1, alpha1, alpha2)

    b_right_left = rbeta(1, beta3, beta4)

    b_right_left_left = rbeta(1, alpha3, alpha4)

    b_right_right_left = rbeta(1, alpha5, alpha6)

    p[i, 1] <- b_left * b_left_left
    p[i, 2] <- b_left * (1 - b_left_left)
    p[i, 3] <- (1 - b_left) * b_right_left * b_right_left_left
    p[i, 4] <- (1 - b_left) * b_right_left * (1 - b_right_left_left)
    p[i, 5] <- (1 - b_left) * (1 - b_right_left) * b_right_right_left
    p[i, 6] <- (1 - b_left) * (1 - b_right_left) * (1 - b_right_right_left)
  }
  return(p)
}


lgnml_approx_mean <- function(beta1, beta2, beta3, beta4,
                              alpha1, alpha2, alpha3, alpha4, alpha5, alpha6){
  gamma1 <- beta1 - alpha1 - alpha2
  gamma2 <- beta3 - alpha3 - alpha4
  gamma3 <- beta4 - alpha5 - alpha6
  gamma4 <- beta2 - beta3 - beta4

  b <- rep(0, 6)

  b[1] <-digamma(alpha1) + digamma(alpha1 + alpha2 + gamma1) - digamma(alpha1 + alpha2) -
    digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  b[2] <- digamma(alpha2) + digamma(alpha1 + alpha2 + gamma1) - digamma(alpha1 + alpha2) -
    digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  b[3] <- digamma(alpha3) + digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + digamma(alpha3 + alpha4 + gamma2) -
    digamma(alpha3 + alpha4) - digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  b[4] <- digamma(alpha4) + digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + digamma(alpha3 + alpha4 + gamma2) -
    digamma(alpha3 + alpha4) - digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  b[5] <- digamma(alpha5) + digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + digamma(alpha5 + alpha6 + gamma3) -
    digamma(alpha5 + alpha6) - digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)


  b[6] <- digamma(alpha6) + digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + digamma(alpha5 + alpha6 + gamma3) -
    digamma(alpha5 + alpha6) - digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    digamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  b <- b - b[6]
  return(b[1:5])
}



lgnml_approx_var <- function(beta1, beta2, beta3, beta4,
                              alpha1, alpha2, alpha3, alpha4, alpha5, alpha6){
  gamma1 <- beta1 - alpha1 - alpha2
  gamma2 <- beta3 - alpha3 - alpha4
  gamma3 <- beta4 - alpha5 - alpha6
  gamma4 <- beta2 - beta3 - beta4


  c_matrix <- matrix(0, ncol = 6, nrow = 6)

  c_matrix[1, 1] <- trigamma(alpha1) + trigamma(alpha1 + alpha2 + gamma1) - trigamma(alpha1 + alpha2) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[2, 2] <- trigamma(alpha2) + trigamma(alpha1 + alpha2 + gamma1) - trigamma(alpha1 + alpha2) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[3, 3] <- trigamma(alpha3) + trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + trigamma(alpha3 + alpha4 + gamma2) -
    trigamma(alpha3 + alpha4) - trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  c_matrix[4, 4] <- trigamma(alpha4) + trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + trigamma(alpha3 + alpha4 + gamma2) -
    trigamma(alpha3 + alpha4) - trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  c_matrix[5, 5] <- trigamma(alpha5) + trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + trigamma(alpha5 + alpha6 + gamma3) -
    trigamma(alpha5 + alpha6) - trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)


  c_matrix[6, 6] <- trigamma(alpha6) + trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + trigamma(alpha5 + alpha6 + gamma3) -
    trigamma(alpha5 + alpha6) - trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)


  ####

  c_matrix[1, 2] <- trigamma(alpha1 + alpha2 + gamma1) - trigamma(alpha1 + alpha2) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[1, 3] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[1, 4] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[1, 5] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[1, 6] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  ##

  c_matrix[2, 3] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[2, 4] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[2, 5] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  c_matrix[2, 6] <- -trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4)

  ##

  c_matrix[3, 4] <- trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + trigamma(alpha3 + alpha4 + gamma2) -
    trigamma(alpha3 + alpha4) - trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  c_matrix[3, 5] <- trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  c_matrix[3, 6] <-trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  ##

  c_matrix[4, 5] <- trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  c_matrix[4, 6] <- trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  ##

  c_matrix[5, 6] <- trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3 + gamma4) + trigamma(alpha5 + alpha6 + gamma3) -
    trigamma(alpha5 + alpha6) - trigamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + gamma1 + gamma2 + gamma3 + gamma4) -
    trigamma(alpha3 + alpha4 + alpha5 + alpha6 + gamma2 + gamma3)

  c_matrix[lower.tri(c_matrix)]  <- t(c_matrix)[lower.tri(c_matrix)]

  ##############

  v_matrix <- matrix(0, ncol = 5, nrow = 5)

  for(i in 1:5){
    for(j in i:5){
      v_matrix[i, j] <- c_matrix[i, j] - c_matrix[i, 6] - c_matrix[j, 6] + c_matrix[6, 6]
    }
  }

  v_matrix[lower.tri(v_matrix)]  <- t(v_matrix)[lower.tri(v_matrix)]

  return(v_matrix)
}



############################################################################################



