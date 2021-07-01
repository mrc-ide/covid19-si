set.seed(42)
library(cowplot)
library(ggforce)
library(ggmcmc)
library(ggplot2)
library(magrittr)
library(matrixStats)
library(patchwork)
source("setup.R")
source("R/utils_process_fits_common.R")
source("R/utils_process_nf_fits.R")
source("R/utils_model_selection.R")
source("R/utils_process_beta_fits.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

nsim_pre_filter <- 20000
nsim_post_filter <- 300
alpha_invalid <- 1
beta_invalid <- 1
max_si <- 40


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

