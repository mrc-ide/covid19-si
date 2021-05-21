set.seed(42)
library(cowplot)
library(dplyr)
library(epitrix)
library(ggforce)
library(glue)
library(ggmcmc)
library(ggplot2)
library(magrittr)
library(matrixStats)
library(patchwork)
library(purrr)
library(rstan)
library(tibble)
source("R/utils.R")
source("R/utils_process_fits_common.R")
source("R/utils_process_nf_fits.R")
source("R/utils_model_selection.R")
source("R/utils_process_beta_fits.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Larger max_shed for NF distribution
max_shed <- 40
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

model_features <- list(
  "mixture" = c(TRUE, FALSE),
  ##"left_bias" = c(TRUE, FALSE),
  "recall"  = c(TRUE, FALSE),
  "right_bias" = c(TRUE, FALSE)
)
model_features <- expand.grid(model_features)
model_features$model_prefix <- ifelse(
  model_features$`right_bias`, "scenario4a", "scenario3a"
)

model_features$model_prefix <-ifelse(
  model_features$mixture,
  glue::glue("{model_features$model_prefix}_mixture"),
  model_features$model_prefix
)

## model_features$model_prefix <-ifelse(
##   model_features$`left_bias`,
##   glue::glue("{model_features$model_prefix}_leftbias"),
##   model_features$model_prefix
## )

model_features$model_prefix <-ifelse(
  model_features$`recall`,
  glue::glue("{model_features$model_prefix}_recall"),
  model_features$model_prefix
)

short_run <- TRUE
iter <- ifelse(short_run, 100, 4000)
chains <- ifelse(short_run, 1, 4)

params_inc <- params_real$inc_par2
si_vec <- seq(-20, max_valid_si)
standata <- list(
  N = nrow(cowling_data), si = cowling_data$si, max_shed = max_shed,
  alpha2 = params_inc[["shape"]], beta2 = 1 / params_inc[["scale"]],
  M = length(si_vec), si_vec = si_vec, width = 0.5
)

#######

fit_model <- function(mixture, recall, right_bias, model_prefix) {
  prefix <- glue("{model_prefix}_nf")
  infile <- glue("stan-models/{prefix}.stan")
  message(infile)
  if (!file.exists(infile)) message("Does not exist ", infile)
  if(mixture) {
    standata$max_invalid_si <- max_invalid_si
    standata$min_invalid_si <- min_invalid_si
  }
  if (right_bias) standata$nu <- cowling_data$nu
  fit <- stan(
    file = infile, data = standata,  verbose = FALSE, iter = iter,
    chains = chains
  )
  outfile <- glue("stanfits/{prefix}_fit.rds")
  saveRDS(fit, outfile)
  fit
}
