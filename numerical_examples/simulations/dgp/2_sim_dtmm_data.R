

sim_dtmm = function(seed,
                    n1 = 40,
                    n2 = 30,
                    n3 = 20,
                    alpha = 1,
                    beta = 0.1,
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

    n_1 = n_2 = n_3= n_4 = n_5 = n_6 = 0

    b_left = rbeta(1, 12 * alpha, 12 * alpha)
    n_left = rbinom(1, N_1[i], b_left)
    n_right = N_1[i] - n_left

    if(n_left > 0){
      b_left_left = rbeta(1, 10 * alpha, 2 * alpha)
      n_1 = rbinom(1, n_left, b_left_left)
      n_2 = n_left - n_1
    }

    if(n_right > 0){
      b_right_left = rbeta(1, 8 * beta, 4 * beta)
      n_right_left = rbinom(1, n_right, b_right_left)
      n_right_right = n_right - n_right_left

      if(n_right_left > 0){
        b_right_left_left = rbeta(1, 4 * beta, 4 * beta)
        n_3 = rbinom(1, n_right_left, b_right_left_left)
        n_4 = n_right_left - n_3
      }
      if(n_right_right > 0){
        b_right_right_left = rbeta(1, 2 * beta, 2 * beta)
        n_5 = rbinom(1, n_right_right, b_right_right_left)
        n_6 = n_right_right - n_5
      }
    }

    y_1[i, ] = c(n_1, n_2, n_3, n_4, n_5, n_6)
  }


  for(i in 1:n2){

    n_1 = n_2 = n_3= n_4 = n_5 = n_6 = 0

    b_left = rbeta(1, 12 * alpha, 12 * alpha)
    n_left = rbinom(1, N_2[i], b_left)
    n_right = N_2[i] - n_left

    if(n_left > 0){
      b_left_left = rbeta(1, 6 * alpha, 6 * alpha)
      n_1 = rbinom(1, n_left, b_left_left)
      n_2 = n_left - n_1
    }

    if(n_right > 0){
      b_right_left = rbeta(1, 8 * beta, 4 * beta)
      n_right_left = rbinom(1, n_right, b_right_left)
      n_right_right = n_right - n_right_left

      if(n_right_left > 0){
        b_right_left_left = rbeta(1, 4 * beta, 4 * beta)
        n_3 = rbinom(1, n_right_left, b_right_left_left)
        n_4 = n_right_left - n_3
      }
      if(n_right_right > 0){
        b_right_right_left = rbeta(1, 2 * beta, 2 * beta)
        n_5 = rbinom(1, n_right_right, b_right_right_left)
        n_6 = n_right_right - n_5
      }
    }

    y_2[i, ] = c(n_1, n_2, n_3, n_4, n_5, n_6)
  }

  for(i in 1:n3){

    n_1 = n_2 = n_3= n_4 = n_5 = n_6 = 0

    b_left = rbeta(1, 12 * alpha, 12 * alpha)
    n_left = rbinom(1, N_3[i], b_left)
    n_right = N_3[i] - n_left

    if(n_left > 0){
      b_left_left = rbeta(1, 2 * alpha, 10 * alpha)
      n_1 = rbinom(1, n_left, b_left_left)
      n_2 = n_left - n_1
    }

    if(n_right > 0){
      b_right_left = rbeta(1, 8 * beta, 4 * beta)
      n_right_left = rbinom(1, n_right, b_right_left)
      n_right_right = n_right - n_right_left

      if(n_right_left > 0){
        b_right_left_left = rbeta(1, 4 * beta, 4 * beta)
        n_3 = rbinom(1, n_right_left, b_right_left_left)
        n_4 = n_right_left - n_3
      }
      if(n_right_right > 0){
        b_right_right_left = rbeta(1, 2 * beta, 2 * beta)
        n_5 = rbinom(1, n_right_right, b_right_right_left)
        n_6 = n_right_right - n_5
      }
    }

    y_3[i, ] = c(n_1, n_2, n_3, n_4, n_5, n_6)
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

