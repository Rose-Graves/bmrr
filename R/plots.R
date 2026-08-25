

#' Plot mu trace and density
#'
#' @param mu_posterior the posterior mu object obtained by running the run_bmrr function
#' @param item_names a vector of items names that you want labeled on the graph
#' @param criteria_names a vector of criteria names that you want labeled on the graph
#'
#' @returns a list of plots
#' @export
#'
#' @examples
#' iter = 1000
#' num_items = 10
#' num_criteria = 3
#'  variation = c(rep(-2,iter), rep(-1, iter), rep(0, iter), rep(1, iter), rep(2, iter))
#' mu_posterior <- array(
#'   rnorm(iter*num_items*num_criteria, mean = 0 + variation, sd = 1),
#'   dim = c(iter,num_items,num_criteria)
#'   )
#' item_names <- paste0("player_",1:num_items)
#' criteria_names <- paste0("time_",1:num_criteria)
#'  mu_plots(mu_posterior, item_names, criteria_names)
mu_plots <- function(mu_posterior,
                     item_names = FALSE,
                     criteria_names = FALSE
){

  #--- Storage for Plots
  plot_list = list()

  #--- Getting number of criteria
  num_criteria = dim(mu_posterior)[3]

  #--- Getting Criteria names sorted
  if(criteria_names[1] != FALSE){
    criteria_values = criteria_names
  }else{
    criteria_values = seq(1, num_criteria, 1)
  }

  ## Trace & Distribution Plots ------------------------------------------------
  # Loop through each criteria
  for(l in 1:num_criteria){


    ### Trace Plot Generation ----------------------------------------------------

    mu_post_trace_df_og = mu_posterior[,,l] %>%
      as.data.frame()

    #--- Add in items names if applicable (Must be in order of the data presented)
    if(item_names[1] == FALSE){
      mu_post_trace_df <- mu_post_trace_df_og%>%
        mutate(
          iter = row_number()
        )%>%
        pivot_longer(
          cols = -iter,
          names_to = "Item",
          values_to = "mu_posterior"
        )%>%
        mutate(
          Item = paste("Item", str_replace(Item, "V", ""))
        )
    }else{
      names(mu_post_trace_df_og) <- item_names
      mu_post_trace_df <- mu_post_trace_df_og%>%
        mutate(
          iter = row_number()
        )%>%
        pivot_longer(
          cols = -iter,
          names_to = "Item",
          values_to = "mu_posterior"
        )
    }

    g_trace <-   ggplot(data = mu_post_trace_df,
                        aes(
                          x = iter,
                          y = mu_posterior,
                          color = Item,
                          group = Item
                        ))+geom_line() +
      ylab(TeX("$\\mu$ Posterior"))+
      xlab("Iteration")+ theme_bw() +
      theme(legend.position="none")


    plot_list = append(plot_list, list(g_trace))


    ### Distribution Plot Generation ------------------------------------------
    mus_post_dist_df_og <- mu_post_trace_df_og

    if(item_names[1] == FALSE){
      names(mus_post_dist_df_og) <- paste("Item", seq(1, dim(mus_post_dist_df_og)[2], 1))
    }else{
      names(mus_post_dist_df_og) <- item_names
    }

    posterior_mean = as.data.frame(round(t(apply(mus_post_dist_df_og, 2, quantile, c(0.025, 0.5, 0.975))), 2))

    posterior_mean = posterior_mean %>%
      tibble::rownames_to_column("Item")%>%
      arrange(`50%`)%>%
      mutate(
        rank = rank(-`50%`)
      )


    mu_graph_data = mus_post_dist_df_og%>%
      pivot_longer(everything(),
                   names_to = "Item",
                   values_to = "Mu Posterior")

    g_dist <- ggplot(data = mu_graph_data%>%
                       group_by(Item)%>%
                       mutate(
                         mean = mean(`Mu Posterior`)
                       )%>%
                       ungroup()%>%
                       arrange(desc(-mean))%>%
                       mutate(Item = fct_reorder(Item, desc(-mean)))
                     ,
                     aes(
                       y = Item,
                       x = `Mu Posterior`,
                       group = Item
                     ))+
      geom_density_ridges(rel_min_height = 0.001, scale = 2)+
      xlab(TeX("Posterior $\\mu")) +
      ylab("Rank Order") +
      theme_bw()


    plot_list = append(plot_list, list(g_dist))


  }

  return(plot_list)

}


#' Plot mu + betaX trace and density
#'
#' @param mu_posterior the posterior mu object obtained by running the run_bmrr function
#' @param beta_posterior the posterior beta object obtained by running the run_bmrr function
#' @param covars the covariate data array
#' @param item_names a vector of items names that you want labeled on the graph
#' @param criteria_names a vector of criteria names that you want labeled on the graph
#'
#' @returns a list of plots
#' @export
#'
#' @examples
#'
#' variation = c(rep(-2,iter), rep(-1, iter), rep(0, iter), rep(1, iter), rep(2, iter))
#' mu_posterior <- array(
#'   rnorm(iter*num_items*num_criteria, mean = 0 + variation, sd = 1),
#'   dim = c(iter,num_items,num_criteria)
#' )
#'
#' beta_posterior <- array(
#'   rnorm(iter*num_covars, mean = 0 + variation, sd = 1),
#'   dim = c(iter,num_covars)
#' )
#'
#' covars <- array(
#'   rnorm(num_items*num_covars*num_criteria, mean = 0, sd = 1),
#'   dim = c(num_items, num_covars, num_criteria)
#' )
#'
#'
#'
#' item_names <- paste("player_",1:num_items)
#' criteria_names <- paste("time_",1:num_criteria)
#' mu_beta_plots(mu_posterior, beta_posterior, covars, item_names, criteria_names)
mu_beta_plots <- function(mu_posterior,
                          beta_posterior,
                          covars,
                          item_names = FALSE,
                          criteria_names = FALSE
){


  #--- Storage for Plots
  plot_list = list()

  #--- Getting number of criteria
  num_criteria = dim(mu_posterior)[3]

  #--- Getting Criteria names sorted
  if(criteria_names[1] != FALSE){
    criteria_values = criteria_names
  }else{
    criteria_values = seq(1, num_criteria, 1)
  }

  ## Trace & Distribution Plots ------------------------------------------------
  # Loop through each criteria
  for(l in 1:num_criteria){


    ### Trace Plot Generation ----------------------------------------------------

    mu_post_trace_df_og = (mu_posterior[,,l] + (beta_posterior%*%t(covars[,,l]))) %>%
      as.data.frame()

    #--- Add in Items names if applicable (Must be in order of the data presented)
    if(item_names[1] == FALSE){
      mu_post_trace_df <- mu_post_trace_df_og%>%
        mutate(
          iter = row_number()
        )%>%
        pivot_longer(
          cols = -iter,
          names_to = "Item",
          values_to = "mu_posterior"
        )%>%
        mutate(
          Item = paste("Item", str_replace(Item, "V", ""))
        )%>%
        filter(
          iter %% 5 == 0
        )
    }else{
      names(mu_post_trace_df_og) <- item_names
      mu_post_trace_df <- mu_post_trace_df_og%>%
        mutate(
          iter = row_number()
        )%>%
        pivot_longer(
          cols = -iter,
          names_to = "Item",
          values_to = "mu_posterior"
        )%>%
        filter(
          iter %% 5 == 0
        )
    }

    g_trace <-   ggplot(data = mu_post_trace_df,
                        aes(
                          x = iter,
                          y = mu_posterior,
                          color = Item,
                          group = Item
                        ))+geom_line() +
      ylab(TeX("$\\mu_{j,k} + \\beta X_{j,k}$ Posterior"))+
      xlab("Iteration")+theme_bw()+
      theme(legend.position="none") +
      xlab(NULL) +
      ylab(NULL)

    plot_list = append(plot_list, list(g_trace))


    ### Distribution Plot Generation ------------------------------------------
    mus_post_dist_df_og <- mu_post_trace_df_og

    if(item_names[1] == FALSE){
      names(mus_post_dist_df_og) <- paste("Item", seq(1, dim(mus_post_dist_df_og)[2], 1))
    }else{
      names(mus_post_dist_df_og) <- item_names
    }

    posterior_mean = as.data.frame(round(t(apply(mus_post_dist_df_og, 2, quantile, c(0.025, 0.5, 0.975))), 2))

    posterior_mean = posterior_mean %>%
      tibble::rownames_to_column("Item")%>%
      arrange(`50%`)%>%
      mutate(
        rank = rank(-`50%`)
      )


    mu_graph_data = mus_post_dist_df_og%>%
      pivot_longer(everything(),
                   names_to = "Item",
                   values_to = "Mu Posterior")

    g_dist <- ggplot(data = mu_graph_data%>%
                       group_by(Item)%>%
                       mutate(
                         mean = mean(`Mu Posterior`)
                       )%>%
                       ungroup()%>%
                       arrange(desc(-mean))%>%
                       mutate(Item = fct_reorder(Item, desc(-mean)))
                     ,
                     aes(
                       y = Item,
                       x = `Mu Posterior`,
                       group = Item
                     ))+
      geom_density_ridges(rel_min_height = 0.001, scale = 2)+
      xlab(TeX("Posterior $\\mu_{j,k} + \\beta X_{j,k}$")) +
      ylab("Rank Order") +
      theme_bw(base_size = 14)


    plot_list = append(plot_list, list(g_dist))


  }

  return(plot_list)

}



#' Plot xi trace and density
#'
#' @param xi_posterior the posterior xi object obtained by running the run_bmrr function
#' @param item_names a vector of items names that you want labeled on the graph
#'
#' @returns a list of plots
#' @export
#'
#' @examples
#' iter = 1000
#' num_items = 10
#' variation = c(rep(-2,iter), rep(-1, iter), rep(0, iter), rep(1, iter), rep(2, iter))
#' xi_posterior <- array(
#'   rnorm(iter*num_items, mean = 0 + variation, sd = 1),
#'   dim = c(iter,num_items)
#'   )
#' item_names <- paste0("player_",1:num_items)
#'  xi_plots(mu_posterior, item_names)
xi_plots <- function(xi_posterior,
                     item_names = FALSE
){

  #--- Storage for Plots
  plot_list = list()


  ## Trace & Distribution Plots ------------------------------------------------

  ### Trace Plot Generation ----------------------------------------------------

  xi_post_trace_df_og = xi_posterior %>%
    as.data.frame()

  #--- Add in Items names if applicable (Must be in order of the data presented)
  if(item_names[1] == FALSE){
    xi_post_trace_df <- xi_post_trace_df_og%>%
      mutate(
        iter = row_number()
      )%>%
      pivot_longer(
        cols = -iter,
        names_to = "Item",
        values_to = "xi_posterior"
      )%>%
      mutate(
        Item = paste("Item", str_replace(Item, "V", ""))
      )
  }else{
    names(xi_post_trace_df_og) <- item_names
    xi_post_trace_df <- xi_post_trace_df_og%>%
      mutate(
        iter = row_number()
      )%>%
      pivot_longer(
        cols = -iter,
        names_to = "Item",
        values_to = "xi_posterior"
      )
  }

  g_trace <-   ggplot(data = xi_post_trace_df,
                      aes(
                        x = iter,
                        y = xi_posterior,
                        color = Item,
                        group = Item
                      ))+geom_line() +
    ylab(TeX("$\\xi$ Posterior"))+
    xlab("Iteration")+
    theme_bw()+
    theme(legend.position="none")

  plot_list = append(plot_list, list(g_trace))


  ### Distribution Plot Generation ------------------------------------------
  xi_post_dist_df_og <- xi_post_trace_df_og

  if(item_names[1] == FALSE){
    names(xi_post_dist_df_og) <- paste("Item", seq(1, dim(xi_post_dist_df_og)[2], 1))
  }else{
    names(xi_post_dist_df_og) <- item_names
  }

  posterior_mean = as.data.frame(round(t(apply(xi_post_dist_df_og, 2, quantile, c(0.025, 0.5, 0.975))), 2))

  posterior_mean = posterior_mean %>%
    tibble::rownames_to_column("Item")%>%
    arrange(`50%`)%>%
    mutate(
      rank = rank(-`50%`)
    )


  xi_graph_data = xi_post_dist_df_og%>%
    pivot_longer(everything(),
                 names_to = "Item",
                 values_to = "Xi Posterior")

  g_dist <- ggplot(data = xi_graph_data%>%
                     group_by(Item)%>%
                     mutate(
                       mean = mean(`Xi Posterior`)
                     )%>%
                     ungroup()%>%
                     arrange(desc(-mean))%>%
                     mutate(Item = fct_reorder(Item, desc(-mean)))
                   ,
                   aes(
                     y = Item,
                     x = `Xi Posterior`,
                     group = Item
                   ))+
    geom_density_ridges(rel_min_height = 0.001, scale = 2)+
    xlab(TeX("Posterior $\\xi")) +
    ylab("Rank Order") +
    theme_bw(base_size = 14)


  plot_list = append(plot_list, list(g_dist))


  return(plot_list)

}


#' Plot beta trace and density
#'
#' @param beta_posterior the posterior beta object obtained by running the run_bmrr function
#' @param covariate_names a vector of covariate names that you want labeled on the graph
#'
#' @returns a list of plots
#' @export
#'
#' @examples
#' iter = 1000
#' num_covars = 5
#' variation = c(rep(-2,iter), rep(-1, iter), rep(0, iter), rep(1, iter), rep(2, iter))
#' beta_posterior <- array(
#'   rnorm(iter*num_covars, mean = 0 + variation, sd = 1),
#'   dim = c(iter,num_covars)
#'   )
#' covariate_names <- paste0("covariate_",1:num_covars)
#' beta_plots(beta_posterior, covariate_names)
beta_plots <- function(beta_posterior,
                       covariate_names = FALSE){


  if(covariate_names[1] != FALSE){
    beta_posterior_og <- data.frame(beta_posterior)
    names(beta_posterior_og) <- covariate_names
  }else{
    beta_posterior_og <- beta_posterior
  }

  beta_posterior_long <- beta_posterior_og %>%
    pivot_longer(
      everything(),
      names_to = "Variable",
      values_to = "BetaPosterior"
    )

  g1 <- ggplot(data = beta_posterior_long,
               aes(x = BetaPosterior,
                   group = Variable,
                   color = Variable
               ))+geom_boxplot()+theme_bw()



  g_dist <- ggplot(data = beta_posterior_long%>%
                     group_by(Variable)%>%
                     filter(Variable != 'Intercept')%>%
                     mutate(
                       mean = mean(`BetaPosterior`)
                     )%>%
                     ungroup()%>%
                     arrange(desc(-mean))%>%
                     mutate(Variable = fct_reorder(Variable, desc(-mean)))
                   ,
                   aes(
                     y = Variable,
                     x = `BetaPosterior`,
                     group = Variable
                   ))+
    geom_density_ridges(rel_min_height = 0.001, scale = 2) + theme_bw(base_size = 16) +
    xlab(TeX("Posterior $\\beta$"))+
    scale_y_discrete(labels = scales::label_parse())

  plot_list = list(g1, g_dist)

  return(plot_list)

}



#' Plot sigma  density
#'
#' @param sigma_posterior the posterior sigma object obtained by running the run_bmrr function
#' @param ranker_names a vector of ranker names that you want labeled on the graph
#'
#' @returns a list of plots
#' @export
#'
#' @examples
#' iter = 1000
#' num_raters = 8
#' variation = c(rep(0,iter), rep(.5, iter), rep(1, iter), rep(2, iter))
#' sigma_posterior <- array(
#'   abs(rnorm(iter*num_raters, mean = 0 + variation, sd = 1)),
#'   dim = c(iter,num_raters)
#'   )
#' ranker_names <- paste0("rater_",1:num_raters)
#' sigma_plots(sigma_posterior, ranker_names)
sigma_plots <- function(sigma_posterior,
                        ranker_names){

  sigma_posterior <- sigma_posterior %>% data.frame()
  names(sigma_posterior) <- ranker_names

  sigma_posterior_long <- sigma_posterior %>%
    pivot_longer(
      everything(),
      names_to = "Rater",
      values_to = "Sigma"
    )

  g1 <- ggplot(data = sigma_posterior_long,
               aes(x = Sigma,
                   group = Rater,
                   color = Rater
               ))+geom_boxplot()


  g_dist <- ggplot(data = sigma_posterior_long%>%
                     group_by(Rater)%>%
                     mutate(
                       `Sigma Posterior` = Sigma,
                       mean = mean(Sigma)
                     )%>%
                     ungroup()%>%
                     arrange(desc(-mean))%>%
                     mutate(Rater = fct_reorder(Rater, desc(-mean)))
                   ,
                   aes(
                     y = Rater,
                     x = `Sigma Posterior`,
                     group = Rater
                   ))+
    geom_density_ridges(rel_min_height = 0.001, scale = 2)+
    xlab(TeX("Posterior $\\sigma^2")) +
    ylab("Raters") +
    theme_bw(base_size = 16)

  #--- Storage for Plots
  plot_list = list(g1,g_dist)

  return(plot_list)
}


#' Plot alpha trace and density
#'
#' @param alpha_posterior the posterior alpha object obtained by running the run_bmrr function
#'
#' @returns a plot
#' @export
#'
#' @examples
#' iter = 1000
#' num_items = 10
#' variation = c(rep(0,iter), rep(.5, iter), rep(1, iter), rep(2, iter), rep(3, iter))
#' alpha_posterior <- array(
#'   abs(rnorm(iter*num_items, mean = 0 + variation, sd = 1)),
#'   dim = c(iter,num_items)
#'   )
#'  alpha_plots(alpha_posterior)
alpha_plots <- function(alpha_posterior){

  alpha_posterior_long <- alpha_posterior %>%
    data.frame()%>%
    pivot_longer(
      everything(),
      names_to = "Item",
      values_to = "Alpha"
    )

  g1 <- ggplot(data = alpha_posterior_long,
               aes(x = Alpha,
                   group = Item,
                   color = Item
               ))+geom_boxplot()+theme_bw()


  return(g1)
}

