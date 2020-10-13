library(dplyr)
library(rstan)
library(hermione)
library(ggmcmc)
library(ggplot2)
library(epitrix)
library(purrr)
library(sbcrs)
source("R/utils.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

max_shed <- 21

## https://doi.org/10.1016/S1473-3099(20)30287-5
mean_inf <- 6.3 ## days
sd_inf <- 4.2 ## days
##params_inf <- epitrix::gamma_mucv2shapescale(mean_inf, sd_inf/ mean_inf)
params_inf <- beta_muvar2shape1shape2(
  mean_inf/max_shed, sd_inf^2 / max_shed^2
)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7014672/
mean_inc <- 6.5 ## days
sd_inc <- 2.6 ## days
params_inc <- epitrix::gamma_mucv2shapescale(mean_inc, sd_inc/ mean_inc)


## original incubation period from IBM
mean_inc_og <- 5.1 #days
sd_inc_og <- 3.94 #days
params_inc_og <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = (sd_inc/mean_inc))

## very short incubation period for testing

# mean_inc <- 1
# sd_inc <- 0.5
# params_inc <- epitrix::gamma_mucv2shapescale(mean_inc, sd_inc/ mean_inc)


mean_iso <- 4
sd_iso <- 2
params_iso <- epitrix::gamma_mucv2shapescale(mean_iso, sd_iso / mean_iso)

