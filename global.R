set.seed(42)
library(dplyr)
library(rstan)
library(hermione)
library(ggmcmc)
library(ggplot2)
library(epitrix)
library(purrr)
library(sbcrs)
library(matrixStats)
library(patchwork)
library(ggforce)
library(tibble)
library(magrittr)
source("R/utils.R")
source("R/scenario_3a_mix_utils.R")
source("R/scenario_4a_mix_utils.R")
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

# additional params for sim data

params <- list(
  inf_par1 = list(mean_inf = 2, sd_inf = 1), ## short
  inf_par2 = list(mean_inf = 8, sd_inf = 2), ## long
  inc_par1 = list(mean_inc = 2, sd_inc = 1), ## short
  inc_par2 = list(mean_inc = 6, sd_inc = 2), ## long
  iso_par1 = list(mean_iso = 2, sd_iso = 2), ## short
  iso_par2 = list(mean_iso = 10, sd_iso = 3), ## long
  offset1 = 0,
  offset2 = -1,
  offset3 = -2,
  pinvalid1 = 0,
  pinvalid2 = 0.05,
  pinvalid3 = 0.2,
  recall1 = 0,
  recall2 = 0.5,
  recall3 = 1
)

# additional parameters for real data

params_real <- list(
  inc_par1 =  list(shape = 5.807, scale = 0.948), # Lauer et al gamma
  inc_par2 = list(shape = 1.675513, scale = 3.043843),# IBM
  inc_par3 = list(shape = 5.807, scale = 0.948), # long - to be changed
  offset1 = -1,
  offset2 = -2,
  offset3 = -3,
  maxshed1 = 9,
  maxshed2 = 14,
  maxshed3 = 21
)


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

## Sanity tests for Neil's distribution implementation
## rstan::expose_stan_functions("stan-models/likelihoods.stan")
## Assume c = 1
## ## a = 0, b = 1, should be 0
## nf_lpdf(0, 1, 1, 10, 5)
## ## a = 1, b = 0, should be 0
## nf_lpdf(1, 0, 1, 10, 5)
## a = b, should be log(1/cosh(a * (t - tmax)))
## nf_lpdf(1, 1, 1, 10, 5)
## 1 / cosh(10 - 5)
