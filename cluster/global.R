set.seed(42)
## library(dplyr)
## library(rstan)
## library(hermione)
## library(ggmcmc)
## library(ggplot2)
## library(epitrix)
## library(purrr)
## library(sbcrs)
## library(matrixStats)
## source("R/utils.R")
## source("R/scenario_3a_mix_utils.R")
##source("R/scenario_4a_mix_utils.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

max_shed <- 21
nsim_pre_filter <- 20000
nsim_post_filter <- 300
alpha_invalid <- 1
beta_invalid <- 1
min_invalid_si <- -11 ## Real has as large as -11
max_si <- 40
max_invalid_si <- 40
max_valid_si <- 40
width <- 0.1


## Further investigation.
##  distribution parameters assume max_shed = 21
params_check <- list(
  inf_par1 = list(mean_inf = 1.5, sd_inf = 0.5), ## short (comparable to s3 real data)
  inf_par2 = list(mean_inf = 5, sd_inf = 5), ## long (comparable to s4 real data)
  inc_par1 = list(shape = 4, scale = 0.5), ## short
  inc_par2 = list(shape = 1.675513, scale = 3.043843),# IBM
  offset1 = -3,
  offset2 = -6,
  pinvalid1 = 0.05,
  pinvalid2 = 0.2,
  beta1 = 0,
  iso_par1 = list(shape = 1, scale = 5)
)
