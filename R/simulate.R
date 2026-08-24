

#' Sample the hierarchical mean values, xi, for each of the items ranked
#'
#' @param num_items The total number of items that were ranked: integer
#'
#' @returns a vector of the hierarchical xi values of length num_items
#' @export
#'
#' @examples
#' get_true_xi(10)
get_true_xi <- function(num_items){
  #--- Xi (hierarchical mean)
  sim_xi = sort(rnorm(num_items, mean = 0, sd = 5))
  return(sim_xi)
}


#' Sample the error values, alpha, for each of the items ranked
#'
#' @param num_items The total number of items that were ranked: integer
#' @param same If you would like the same error value for each item TRUE else FALSE
#'
#' @returns a vector of alpha values of length num_items
#' @export
#'
#' @examples
#' get_true_alpha(10, TRUE)
get_true_alpha <- function(num_items, same){
  #--- Alpha (hierarchical variance)
  if(same == TRUE){
    sim_alpha <- rep(.1, num_items)
    return(sim_alpha)
  }else{
    sim_alpha <- rgamma(num_items, shape = 2, scale = 1/2)
    return(sim_alpha)
  }
}


#' Sample the mean values for each of the items for each criteria
#'
#' @param xi_true a vector of length num_items that is the true xi values
#' @param xi_alpha a vector of length num_items that is the true alpha values
#' @param num_items the total number of items ranked: integer
#' @param num_criteria the number of criteria/time points: integer
#'
#' @returns a matrix of mu values for each item across each criteria/time point
#' @export
#'
#' @examples
#' num_items = 10
#' num_criteria = 3
#' xi_true <- get_true_xi(num_items)
#' xi_alpha <- get_true_alpha(num_items, TRUE)
#' get_true_mu(xi_true, xi_alpha, num_items, num_criteria)
get_true_mu <- function(xi_true, xi_alpha, num_items, num_criteria){
  #--- Mu (random effects)
  sim_mu <- matrix(
    data = NA,
    nrow = num_items,
    ncol = num_criteria
  )
  for(j in 1:num_items){
    sim_mu_j <- rnorm(num_criteria, xi_true[j], xi_alpha[j])
    sim_mu[j,] <- sim_mu_j
  }
  return(sim_mu)
}


#' Sample the rater error
#'
#' @param num_raters
#' @param same If you would like the same error value for each rater TRUE else FALSE
#' @param sigma_value If you want the same error value, specify it here
#'
#' @returns
#' @export
#'
#' @examples
#' get_sigma_true(8, TRUE, 1)
get_sigma_true <- function(num_raters, same, sigma_value = 1){
  #--- epsilon (rater error)
  if(same == TRUE){
    sim_sigma <- rep(sigma_value, num_raters)
    return(sim_sigma)
  }else{
    sim_sigma <- rgamma(num_raters, shape = 2, scale = 1/2)
    return(sim_sigma)
  }
}


#' Sample the coefficient values
#'
#' @param num_covars the number of covariates: integer
#'
#' @returns a vector of beta values
#' @export
#'
#' @examples
#' get_true_beta(5)
get_true_betas <- function(num_covars){
  #--- Betas
  sim_betas <- runif(num_covars,-5,5)
  return(sim_betas)
}


#' Simulate the covariates
#'
#' @param num_items the total number of items ranked: integer
#' @param num_covars the number of criteria/time points: integer
#'
#' @returns
#' @export
#'
#' @examples
#' simulate_covars(10, 5)
simulate_covars <- function(num_items, num_covars){
  #--- Covariates
  intercept = rep(1, num_items)
  covar_moment = data.frame(intercept)
  for(c in 1:num_covars){

    covar_now = runif(num_items, -5, 5)

    #exclude minus one to one
    covar_moment <- cbind(covar_moment,  covar_now)
  }
  covar_sim <- repair_names(covar_moment)
  return(covar_sim)
}


#' Simulate multivariate latent data
#'
#' @param player_random_matrix matrix of the latent values for each item over each time point
#' @param sigma_true the values of the rater error
#' @param num_items the total number of items ranked: integer
#' @param num_raters the number of raters: integer
#' @param num_criteria the number of criteria/time points: integer
#' @param partial_raters the number of raters with partial rankings: integer
#' @param missing_percent the percentage of items that partial raters did not rank: real number between 0,1
#' @param rank_R
#'
#' @returns
#' @export
#'
#' @examples
#'
#'
#' num_items = 10
#' num_raters = 8
#' num_criteria = 3
#' partial_raters = 4
#' missing_percent = 0.20
#' sigma_true = get_sigma_true(num_raters, TRUE, 1)
#' player_random_matrix <- matrix( data = 1, nrow = num_items, ncol = num_criteria)
#' simulate_z(player_random_matrix,sigma_true, num_items, num_raters, num_criteria, partial_raters, missing_percent, rank_R = TRUE )
simulate_z <- function(player_random_matrix, sigma_true,
                       num_items, num_raters, num_criteria,
                       partial_raters, missing_percent,
                       rank_R = TRUE){

  #--- define array
  result <- array(NA,dim = c(num_raters, num_items, num_criteria))

  #--- Calculating number of partial rankings
  partial_number <- round(num_raters*missing_percent)
  raters_partial <- sample(1:num_raters, partial_raters, replace= FALSE)

  # simulate data
  for (i in 1:num_raters){

    latent_Z_matrix = matrix(data = NA, nrow = num_items, ncol = num_criteria)

    for(j in 1:num_items){

      for(k in 1:num_criteria){
        # get latent continuous Z
        latent_Z_matrix[j,k] = player_random_matrix[j,k] + rnorm(1, 0 , sigma_true[i])

      }

    }

    if(rank_R){
      for(k in 1:num_criteria){
        result[i,,k] = rank(-latent_Z_matrix[,k])

        if(i %in% raters_partial){
          columns_remove<- which(rank(-latent_Z_matrix[,k]) %in% (num_raters-partial_number+1):num_raters)
          result[i,columns_remove,k] <- NA
        }
      }
    }
    else{ ### Update needed if want this ####
      # else record the latent_Z
      result[i,] = latent_Z
    }
  }
  return(result)


}




#' Simulate multivariate ranked data
#'
#' @param num_items the total number of items ranked: integer
#' @param same_alpha If you would like the same error value for each item TRUE else FALSE
#' @param num_criteria the number of criteria/time points: integer
#' @param num_raters the number of raters: integer
#' @param same_sigma If you would like the same error value for each rater TRUE else FALSE
#' @param num_covars the number of covariates: integer
#' @param partial_raters the number of raters with partial rankings: integer
#' @param missing_percent the percentage of items that partial raters did not rank: real number between 0,1
#' @param seed the seed duh
#' @param sigma_value the amount of error
#'
#' @returns
#' @export
#'
#' @examples
#'  num_items = 10
#'  num_raters = 8
#'  num_criteria =  3
#'  num_covars = 2
#'  covars_present = TRUE
#'  same_alpha = TRUE
#'  same_sigma = TRUE
#'  partial_raters = 4
#'  missing_percent = 0.3
#'  simulate_data <- function(num_items, same_alpha, num_criteria, num_raters,
#'  same_sigma, num_covars, partial_raters, missing_percent, seed, sigma_value = 1)
simulate_data <- function(num_items, same_alpha, num_criteria, num_raters,
                         same_sigma, num_covars, partial_raters,
                         missing_percent, seed, sigma_value = 1){


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


