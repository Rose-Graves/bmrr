
## Overview

The **bmrr** package fits a Bayesian Multivariate Rank Regression
(BMRR). The BMRR model serves to aggregate multivariate ranking data in
cases where raters rank multiple items over multiple time points or
criteria. Using latent hierarchical regression models, BMRR pools
information across raters and criteria, incorporates covariates and
rater expertise, and constructs an aggregate ranking. Crucially, BMRR
provides full posterior uncertainty quantification for this aggregate
ranking and all model parameters, and readily handles partial (i.e.,
incomplete) rankings and ties. In addition to fitting the model, the
package provides helpful posterior distribution plots to analyze aspects
of the aggregated ranking, the raters themselves, the covariate effects,
and provides comparisons across items.

## How to use bmrr

``` r
#--- Load in package
library(bmrr)

#--- Load in some additional packages
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
library(tidyr)
```

Below we walk through an applied example with Major League Baseball
(MLB) ranking data.

``` r
#--- Load in the data
data(sports_data)

# Define Simulation Starting values & Priors ----------------------------------

# Sampling Specifications ----------------------------------------------------
#--- Set seed
seed = 24
set.seed(seed)

#--- Number of samples
samps = 20000

#--- Burn-in
burnin = 4000

#--- number of items 
num_items = m = length(unique(sports_data$player_id))

#--- number of raters
num_raters = n = 5

#--- number of times
num_criteria = t = length(unique(sports_data$Year))

#--- number of covariates
c = 3

#--- Priors for beta and lambda
lambda = beta = .8

#--- Yes Covars are Present
covars_present = TRUE


# Clean Data & Parameters --------------------------------------------

covars = array(NA,dim = c(m,c,t))

for(i in 1:t){
  covars[,,i] <- sports_data %>%
    arrange(player_id)%>%
    filter(Year == 2020+i)%>%
    ungroup()%>%
    dplyr::select(
      player_age, 
      WAR, 
      base_salary,
    )%>%as.matrix()
}

sports_data_array <- array(NA,dim = c(n,m,t))
for(i in 1:t){
  sports_data_array[,,i] <- sports_data%>%
    arrange(player_id)%>%
    filter(Year == 2020+i)%>%
    ungroup()%>%
    dplyr::select(
      ESPN, 
      CBS,
      MLB, 
      BR, 
      Yahoo
    )%>%t()
}

data <- sports_data_array
```

``` r
# Priors -----------------------------------------------------------------------

#--- Prior for xi (hierarchical mean)
mu_0 = rep(0, m) #0
sigma_0 = rep(1,m) #1

#--- Prior for alpha (hierarchical variance)
alpha_0 <- .8

#--- Prior for beta
beta_0 = rep(0, c) #there are 'c' coefficients 
sum_0 = diag(100,c) #there are 'c' coefficients
sum_0_inverse = solve(sum_0)

#--- Prior for Sigma
lambda = .8


# Run Simulation ------------------------------------------------------------
posterior_output <- run_bmrr(num_raters,
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
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%[1] "It's done honeybun!"


# Simulation Results------------------------------------------------------------

#--- Mu Posterior
mu_posterior <- posterior_output$mu_posterior

#--- Xi Posterior
xi_posterior <- posterior_output$xi_posterior

#--- Sigma Posterior
sigma_posterior <- posterior_output$sigma_posterior

#--- Alpha Posterior
alpha_posterior <- posterior_output$alpha_posterior 

#--- Beta Posterior
beta_posterior <- posterior_output$beta_posterior

#--- Z Posterior
Z_posterior <- posterior_output$Z_posterior

#--- Run Time
time <- as.numeric(posterior_output$time)
```

``` r

mu_beta_plot_list <- mu_beta_plots(mu_posterior, 
                         beta_posterior,
                         covars,
                         item_names, 
                         criteria_names)

mu_beta_plot_list[[5]]
mu_beta_plot_list[[6]] # MAIN PLOT FOR PAPER
```

``` r

## Posterior Xi ----------------------------------------------------------------
xi_plot_list <- xi_plots(xi_posterior, 
                         item_names)

xi_plot_list[[1]] 
xi_plot_list[[2]] # PLOT FOR PAPER
```

``` r

covariate_names <- c( "Age", "WAR", "Salary")


beta_plots_list <- beta_plots(beta_posterior, 
                              covariate_names)

beta_plots_list[[2]] # plot for paper
```

``` r
ranker_names <-c("ESPN","CBS", "MLB","BR","Yahoo")
colMeans(sigma_posterior)
colSds(sigma_posterior)

sigma_plots_list <- sigma_plots(sigma_posterior, 
                                ranker_names)
sigma_plots_list
```

``` r

time1 = (mu_posterior[, ,4] + (beta_posterior%*%t(covars[,,4]))) %>% 
  as.data.frame()%>%
  mutate(
    iter = row_number()
  )

compare_prob_pos_table = matrix(NA, nrow = m, ncol = m)
#--- Reading Column to Row (Probability Column 1 is larger than Row 2)

for(j in 1:m){
  for(k in 1:m){
    prob_now = sum(time1[,j] >= time1[,k])/dim(time1)[1]
    if(j == k){
      compare_prob_pos_table[k,j] = 1
    }
    else{
      compare_prob_pos_table[k,j]<- prob_now
    }
  }
}


compare_prob_pos_table <- data.frame(compare_prob_pos_table)%>%
  mutate_if(is.character,as.numeric)

names(compare_prob_pos_table) <-item_names

rownames(compare_prob_pos_table) <- item_names

compare_prob_pos_table = round(compare_prob_pos_table, 3)

compare_prob_pos_table_OF = compare_prob_pos_table %>% 
  dplyr::select(

    #--- Outfield
    "George Springer",
    "Mike Trout",          
    "Bryce Harper",
    "Yordan Alvarez",
    "Aaron Judge",
    "Mookie Betts",
    "Kyle Tucker",
    "Christian Yelich",
    "Teoscar Hernandez" ,
    "Brandon Nimmo",
    "Ronald Acuna Jr",
    "Fernando Tatis Jr",
    "Juan Soto" ,
    "Randy Arozarena" ,
    "Luis Robert Jr" 

  )%>%
  filter(
    row.names(compare_prob_pos_table) %in% c("George Springer",
                       "Mike Trout",          
                       "Bryce Harper",
                       "Yordan Alvarez",
                       "Aaron Judge",
                       "Mookie Betts",
                       "Kyle Tucker",
                       "Christian Yelich",
                       "Teoscar Hernandez" ,
                       "Brandon Nimmo",
                       "Ronald Acuna Jr",
                       "Fernando Tatis Jr",
                       "Juan Soto" ,
                       "Randy Arozarena" ,
                       "Luis Robert Jr" )
  )


compare_prob_pos_table_OF <- (as.matrix(compare_prob_pos_table_OF))


sorted_row_names <- sort(rownames(compare_prob_pos_table_OF))
sorted_col_names <- rev(sort(colnames(compare_prob_pos_table_OF)))


long_data <- melt(compare_prob_pos_table_OF, varnames = c("RowItem", "ColumnItem"), value.name = "Probability")

long_data$RowItem <- factor(long_data$RowItem, levels = sorted_row_names)
long_data$ColumnItem <- factor(long_data$ColumnItem, levels = sorted_col_names)

#--- Create the Heatmap

long_data = long_data %>% 
  mutate(Probability = 1- Probability)%>%
  mutate(Probability = ifelse(RowItem == ColumnItem, 1, Probability))

# Create the base plot
heatmap_plot <- ggplot(
  data = long_data,
  aes(x = ColumnItem, y = RowItem, fill = Probability)
) +
  geom_tile(color = "white")+  
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )+ 
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = .5,
    limit = c(0, 1)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

heatmap_plot
```

## Installation

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("Rose-Graves/bmrr")
```
