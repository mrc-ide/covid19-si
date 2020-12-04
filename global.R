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
source("R/utils.R")
source("R/scenario_1a_mix_utils.R")
source("R/scenario_2a_mix_utils.R")
source("R/scenario_3a_mix_utils.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

max_shed <- 21
nsim <- 500
alpha_invalid <- 1
beta_invalid <- 1
min_si <- -2 ## For simulations
offset <- -2 # this is not the smallest possible SI, but the maximum
## time before symptom onset when secondary infection can happeb.

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
  mu = mean_inc_og, cv = (sd_inc_og/mean_inc_og))

## very short incubation period for testing

# mean_inc <- 1
# sd_inc <- 0.5
# params_inc <- epitrix::gamma_mucv2shapescale(mean_inc, sd_inc/ mean_inc)


mean_iso <- 4
sd_iso <- 2
params_iso <- epitrix::gamma_mucv2shapescale(mean_iso, sd_iso / mean_iso)

params <- list(
  inf_par1 = list(mean_inf = 7.14, sd_inf = 2.94),
  inf_par2 = list(mean_inf = 4.2, sd_inf = 2.73),
  inf_par3 = list(mean_inf = 3.15, sd_inf = 2.31),
  inc_par1 = list(mean_inc = 6.5, sd_inc = 2.6),
  inc_par2 = list(mean_inc = 5.1, sd_inc = 3.94),
  ## very short incubation period to allow simulation of -ve SIs
  inc_par3 = list(mean_inc = 2, sd_inc = 1),
  iso_par1 = list(mean_iso = 2, sd_iso = 2),
  iso_par2 = list(mean_iso = 4, sd_iso = 2),
  iso_par3 = list(mean_iso = 4, sd_iso = 4),
  offset1 = -1,
  offset2 = -2,
  offset3 = -3,
  pinvalid1 = 0.01,
  pinvalid2 = 0.03,
  pinvalid3 = 0.05,
  recall1 = 0.01,
  recall2 = 0.1,
  recall3 = 0.1
)

## Set up params grid
## param_grid <- expand.grid(
##   params_inf = c("inf_par1", "inf_par2", "inf_par3"),
##   params_inc = c("inc_par1", "inc_par2"),
##   params_pinv = c("pinvalid1", "pinvalid2", "pinvalid3"),
##   recall = c("recall1", "recall2", "recall3"),
##   stringsAsFactors = FALSE
## )


