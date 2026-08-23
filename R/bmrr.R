
# Function: Simulating Data ----------------------------------------------------

simulate_data <- function(num_items,
                          same_alpha,
                          num_criteria,
                          num_raters,
                          same_sigma,
                          num_covars,
                          partial_raters,
                          missing_percent,
                          seed,
                          sigma_value = 1
){

  #num_items: The number of items (integer)
  #same_alpha: Do you want the alpha parameter across all items (boolean)
  #num_criteria: The number of criteria (integer)
  #num_raters: The number of raters (integer)
  #same_sigma: Do you want the same sigma parameter across all items (boolean)
  #num_covars: The number of covariates (integer)
  #partial_raters: the number of raters with partial rankings (integer)
  #missing_percent: the percentage of items that partial raters did not rank (real number between 0,1)
  #seed: the seed duh

  #--- Set seed
  set.seed(seed)

  #--- get the hierarchical means
  xi_true = get_true_xi(num_items)

  #--- get hierarchical variance
  alpha_true = get_true_alpha(num_items, same_alpha)

  #--- get random effects
  mu_true = get_true_mu(xi_true, alpha_true, num_items, num_criteria)

  #--- get the rater variance
  sigma_true = get_sigma_true(num_raters, same_sigma, sigma_value)


  if(num_covars >0){

    ##---  Covariates
    #--- Get true betas
    beta_true <- c(1, get_true_betas(num_covars))

    #--- simulate the covariates
    covars <- simulate_covars(num_items,num_covars)

    #--- putting covariates into easier format to use
    X = covars
    for(i in 1:(num_criteria-1)){
      X = X %>% rbind(covars)
    }
    xtx = as.matrix(t(X))%*%as.matrix(X)

    #--- covariate impact on latent z
    covars_impact <- t(beta_true%*%t(as.matrix(covars)))

    #--- replicating covariates for all criteria
    covars_impact_matrix <- matrix(rep(covars_impact, num_criteria), ncol = num_criteria)

    #--- Getting latent z before error
    mu_covars_true = mu_true + covars_impact_matrix

    player_random_matrix = mu_covars_true

    ##--- Simulate the data
    data <- simulate_z(
      player_random_matrix, sigma_true,
      num_items, num_raters, num_criteria,
      partial_raters, missing_percent,
      rank_R = TRUE
    )

    ##--- Organize Data and parameters into list
    sim_data_params <-  list(
      data = data,
      covars = covars,
      xi_true = xi_true,
      alpha_true = alpha_true,
      mu_true = mu_true,
      sigma_true = sigma_true,
      beta_true = beta_true
    )

  }else{

    player_random_matrix = mu_true

    ##--- Simulate the data
    data <- simulate_z(
      player_random_matrix, sigma_true,
      num_items, num_raters, num_criteria,
      partial_raters, missing_percent,
      rank_R = TRUE
    )

    covars = FALSE
    beta_true = FALSE

    ##--- Organize Data and parameters into list
    sim_data_params <-  list(
      data = data,
      covars = covars,
      xi_true = xi_true,
      alpha_true = alpha_true,
      mu_true = mu_true,
      sigma_true = sigma_true,
      beta_true = beta_true
    )

  }


  return(
    sim_data_params
  )

}



# Function: Running Model ------------------------------------------------------

run_bmmr <- function(num_raters,
                     num_items,
                     num_criteria,
                     data,
                     samps,
                     burnin,
                     covars_present,
                     covars,

                     mu_0,
                     sigma_0,
                     alpha_0,
                     beta_0,
                     sum_0_inverse,
                     lambda,
                     seed
){


  #num_raters: number of raters (integer)
  #num_items: number of items (integer)
  #num_criteria: number of criteria (integer)
  #data: The data (array)
  #samps: number of samples post-burin (integer)
  #burnin: number of burnin samples (integer)
  #covars: covariate data

  #--- Prior for xi (hierarchical mean)
  # mu_0: prior mean
  # sigma_0: prior variance

  #--- Prior for alpha (hierarchical variance)
  #alpha_0

  #--- Prior for beta
  #beta_0: prior mean
  #sum_0_inverse: prior variance

  #--- Prior for Sigma
  #lambda

  set.seed(seed)

  n = num_raters
  m = num_items
  K = num_criteria

  #--- Defining intercept for no-covariate situations
  if(covars_present == FALSE){
    covars = array(rep(1, num_items),dim = c(num_items,1,num_criteria))
    X = t(rbind(covars))
    xtx = as.matrix(t(X))%*%as.matrix(X)
  }else{
    #--- Getting X into a long matrix
    X = do.call(rbind, lapply(1:num_criteria, function(f) covars[,,f]))
    xtx = as.matrix(t(X))%*%as.matrix(X)
  }

  #--- Getting number of covariates
  num_covars = dim(covars)[2] - 1

  #--- Initial Z values
  z_ij_last = array(NA, dim = c(num_raters,num_items,num_criteria))

  for(i in 1:num_raters){
    for(l in 1:num_criteria){
      values = -sort(runif(m, -10, 10))
      z_ij_last[i,,l] <- values[data[i,,l]]
    }
  }
  #---- Accounting for any unranked items
  z_ij_last[is.na(z_ij_last)] <- -10

  #--- Getting the totality of Z
  Z_posterior = array(NA, dim = c(samps+burnin, n*m, num_criteria))

  #--- Mu Posterior Values & Initial Values
  mu_posterior = array(NA, dim = c(samps+burnin+1,m,num_criteria))
  sort_indices <- rank(-colMeans(data[,,1]))
  mu_posterior[1,,] = sort(runif(m, -10, 10))[sort_indices]


  #--- Sigma Posterior Values & Initial Values
  sigma_posterior = matrix(1, ncol = n, nrow = samps+burnin+1)
  for(i in 1:(samps+burnin+1)){
    sigma_posterior[i,] = rep(1,n)
  }

  #--- Zeta Posterior Values & Initial Values
  xi_posterior = matrix(NA, ncol = m, nrow = samps+burnin+1)
  xi_posterior[1,] <- sort(runif(m, -10, 10))[sort_indices]

  #--- Alpha Posterior Values & Initial Values
  alpha_posterior = matrix(NA, ncol = m, nrow = samps+burnin+1)
  alpha_posterior[1,] <- rep(1/8, m)

  #--- Setting beta_posterior matrix and setting initial values at 0
  beta_posterior = matrix(NA, ncol = num_covars+1, nrow = samps+burnin+1)
  beta_posterior[1,] = rep(0, num_covars+1)

  pb <- txtProgressBar(min = 1, max = samps+burnin, style = 3)

  start_time <- Sys.time()
  for(k in 1:(samps+burnin)){

    setTxtProgressBar(pb, k)

    ### Step 1: Sampling for values of z_ij -----------------------------------
    for(i in 1:n){
      for(l in 1:K){
        #--- Get ranks from rater i
        current_ranks = data[i,,l]

        #--- determine if partial rank
        is_partial = any(is.na(current_ranks))
        unranked_loc = which(is.na(current_ranks))

        if(is_partial == FALSE){
          for(j in 1:m){
            current_rank_loc = which(current_ranks==j) #Getting placement of item
            if(j==1){
              # Item Ranked First
              lower_rank = which(current_ranks==j+1) #Getting placement item ranked 2nd
              current_mean = mu_posterior[k,current_rank_loc,l]+ sum(covars[current_rank_loc,,l]*beta_posterior[k,]) # Getting the mean

              z_ij = rtruncnorm(
                1,
                a = z_ij_last[i,lower_rank,l],
                b = 10,
                mean = current_mean,
                sd = sqrt(sigma_posterior[k,i])

              )
              z_ij_last[i,current_rank_loc,l] = z_ij

            }
            if(j==m){
              # Item Ranked Last
              upper_rank = which(current_ranks== j-1)
              current_mean = mu_posterior[k,current_rank_loc,l]+ sum(covars[current_rank_loc,,l]*beta_posterior[k,]) # Getting the mean

              z_ij = rtruncnorm(
                1,
                a = -10,
                b = z_ij_last[i,upper_rank,l],
                mean = current_mean,
                sd = sqrt(sigma_posterior[k,i])
              )
              z_ij_last[i,current_rank_loc,l] = z_ij

            }
            if(j!=1 & j!=m){
              # Any Other Item
              lower_rank = which(current_ranks==j+1)
              upper_rank = which(current_ranks==j-1)
              current_mean = mu_posterior[k,current_rank_loc,l]+ sum(covars[current_rank_loc,,l]*beta_posterior[k,]) # Getting the mean

              z_ij = rtruncnorm(
                1,
                a = z_ij_last[i,lower_rank,l],
                b = z_ij_last[i,upper_rank,l],
                mean = current_mean,
                sd = sqrt(sigma_posterior[k,i])
              )
              z_ij_last[i,current_rank_loc,l] = z_ij

            }
          }
        }
        #--- Function for partial ranks
        if(is_partial == TRUE){
          max_rank = max(current_ranks, na.rm = TRUE)
          for(j in 1:max_rank){
            current_rank = which(current_ranks==j)
            if(j==1){
              lower_rank = which(current_ranks==j+1)
              current_mean = mu_posterior[k,current_rank,l]+ sum(covars[current_rank,,l]*beta_posterior[k,]) # Getting the mean

              z_ij = rtruncnorm(
                1,
                a = z_ij_last[i,lower_rank,l],
                b = 10,
                mean = current_mean,
                sd = sqrt(sigma_posterior[k,i])
              )
              z_ij_last[i,current_rank,l] = z_ij

            }
            if(j==max_rank){
              upper_rank = which(current_ranks== j-1)
              current_mean = mu_posterior[k,current_rank,l]+ sum(covars[current_rank,,l]*beta_posterior[k,]) # Getting the mean

              z_ij = rtruncnorm(
                1,
                a = ifelse(is.na(max(z_ij_last[i,unranked_loc,l] )),-10, max(z_ij_last[i,unranked_loc,l] )),
                b = z_ij_last[i,upper_rank,l],
                mean = current_mean,
                sd = sqrt(sigma_posterior[k,i])
              )
              z_ij_last[i,current_rank,l] = z_ij

            }
            if(j!=1 & j<max_rank){
              lower_rank = which(current_ranks==j+1)
              upper_rank = which(current_ranks==j-1)
              current_mean = mu_posterior[k,current_rank,l]+ sum(covars[current_rank,,l]*beta_posterior[k,]) # Getting the mean

              z_ij = rtruncnorm(
                1,
                a = z_ij_last[i,lower_rank,l],
                b = z_ij_last[i,upper_rank,l],
                mean = current_mean,
                sd = sqrt(sigma_posterior[k,i])
              )
              z_ij_last[i,current_rank,l] = z_ij

            }
          }
          for(j in 1:length(unranked_loc)){
            current_rank = unranked_loc[j]
            upper_rank = which(current_ranks==max_rank)
            current_mean = mu_posterior[k,current_rank,l]+ sum(covars[current_rank,,l]*beta_posterior[k,]) # Getting the mean

            z_ij = rtruncnorm(
              1,
              a = -10,
              b = z_ij_last[i,upper_rank,l],
              mean = current_mean,
              sd = sqrt(sigma_posterior[k,i])
            )
            z_ij_last[i,current_rank,l] = z_ij

          }

        }
      }
    }

    #--- Save off standardized values of Z
    for(l in 1:K){
      z_ij_last[,,l] = (z_ij_last[,,l] - mean(z_ij_last[,,l]))/sd(z_ij_last[,,l])
      Z_posterior[k,,l] = c(z_ij_last[,,l])
    }


    ### Step 2: Update Xi ------------------------------------------------------
    xi_samp = c()
    for(j in 1:m){
      mu_j = c()
      for(l in 1:num_criteria){ # Getting all mu_j for each criteria
        mu_j = c(mu_j, mu_posterior[k,j,l])
      }
      xi_mean_numerator = (sum(mu_j)/alpha_posterior[k,j]) + mu_0/sigma_0
      xi_mean_denomerator = (num_criteria/alpha_posterior[k,j]) + 1/sigma_0
      xi_var  = 1/xi_mean_denomerator
      xi_samp = c(xi_samp, rnorm(1, mean = xi_mean_numerator*xi_var, sd = sqrt(xi_var)))

    }
    xi_posterior[k+1,] <- xi_samp

    ### Step 3: Update Alpha --------------------------------------------------
    alpha_samp = c()
    for(j in 1:m){
      mu_zeta_diff = c()
      for(l in 1:K){
        mu_zeta_diff = c(mu_zeta_diff, (mu_posterior[k,j,l]- xi_posterior[k+1,j])^2)
      }
      alpha_alpha = (K/2) + (1/alpha_0)
      alpha_beta = (1/2)*sum(mu_zeta_diff) + alpha_0
      alpha_samp = c(alpha_samp, rinvgamma(1, alpha_alpha, alpha_beta))
    }

    alpha_posterior[k+1,] <- alpha_samp



    ### Step 4: Update Mu --------------------------------------------------
    #--- calculating common values
    sigma_sum_list = c()
    for(i in 1:n){
      sigma_sum_list = c(sigma_sum_list, 1/sigma_posterior[k,i])
    }
    sigma_sum = sum(sigma_sum_list)

    #--- getting posterior of latent mu values for each item
    for(j in 1:m){
      sigma_mu_post = sigma_sum + 1/alpha_posterior[k+1,j]

      for(l in 1:K){
        z_sum_list = c()

        cov_impact_jl = sum(covars[j,,l] * beta_posterior[k,])

        for(i in 1:n){
          z_sum_list = c(z_sum_list, (z_ij_last[i,j,l]- cov_impact_jl) / sigma_posterior[k,i])
        }
        z_sum = sum(z_sum_list, xi_posterior[k+1,j]/alpha_posterior[k+1,j])
        var_mu_post = 1/sigma_mu_post
        mu_jt = rnorm(1, mean = z_sum*var_mu_post, sd = sqrt(var_mu_post))
        mu_posterior[k+1,j,l] = mu_jt
      }
    }


    ### Step 5: Update Beta ---------------------------------------------------
    if(covars_present == TRUE){
      #--- getting sigma & mean portion
      sigma_sum = matrix(0, nrow = dim(xtx)[1], ncol = dim(xtx)[2])
      mean_sum = matrix(0, nrow = dim(xtx)[1], ncol = 1)
      mu_sum = matrix(0, nrow = dim(xtx)[1], ncol = 1)

      for(i in 1:n){
        sigma_sum = sigma_sum + xtx/sigma_posterior[k,i]
        mean_sum = mean_sum + t(X)%*%c(z_ij_last[i,,])/sigma_posterior[k,i]
        mu_sum = mu_sum + (t(X)%*%c(mu_posterior[k+1,,]))/sigma_posterior[k,i]
      }
      sigma_betas = sigma_sum + sum_0_inverse
      sigma_betas_inv = solve(sigma_betas)
      mean_betas =  (mean_sum - mu_sum) + sum_0_inverse%*%beta_0

      #--- sampling portion
      betas_current <- mvrnorm(n = 1, mu = sigma_betas_inv%*%mean_betas, Sigma = sigma_betas_inv)

      beta_posterior[k+1,] <- betas_current
    }else{
      betas_current <- 0
      beta_posterior[k+1,] <- betas_current
    }

    ### Step 6: Gibbs for Sigma --------------------------------------------------
    alpha_sigma_post = (K*m)/2 + 1/lambda

    for(i in 1:n){
      beta_sigma_list = c()

      for(j in 1:m){
        for(l in 1:K){

          sum_squares_value_now = (z_ij_last[i,j,l] - mu_posterior[k+1, j,l] - beta_posterior[k+1,]%*%covars[j,,l]  )^2
          beta_sigma_list = c(beta_sigma_list, sum_squares_value_now)
        }
      }
      beta_sigma_sum = sum(beta_sigma_list)
      beta_sigma_post = (1/2)*beta_sigma_sum + lambda

      #--- Sample Posterior sigma
      sigma_i = rinvgamma(1, alpha_sigma_post, beta_sigma_post)

      #--- save off posterior sigma
      sigma_posterior[k+1, i] = sigma_i
    }



  }
  end_time <- Sys.time()
  time = end_time - start_time


  model_posterior <- list(
    mu_posterior = mu_posterior[(burnin+1):(burnin+samps),,],
    xi_posterior = xi_posterior[(burnin+1):(burnin+samps),],
    sigma_posterior = sigma_posterior[(burnin+1):(burnin+samps),],
    alpha_posterior = alpha_posterior[(burnin+1):(burnin+samps),],
    beta_posterior = beta_posterior[(burnin+1):(burnin+samps),],
    Z_posterior = Z_posterior[(burnin+1):(burnin+samps),,],
    time = time
  )

  print("It's done honeybun!")

  return(model_posterior)

}

