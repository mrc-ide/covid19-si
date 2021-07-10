library(dplyr)
library(epitrix)
library(glue)
library(purrr)
library(rstan)
library(tibble)
source("cowling-data-prep.R")
source("R/utils.R")
max_shed <- 21
min_valid_si <- -20
max_valid_si <- 40
min_invalid_si <- -20
max_invalid_si <- 40
width <- 0.5

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


model_features <- list(
  "mixture" = c(TRUE, FALSE),
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

model_features$model_prefix <-ifelse(
  model_features$`recall`,
  glue::glue("{model_features$model_prefix}_recall"),
  model_features$model_prefix
)

short_run <- FALSE
iter <- ifelse(short_run, 100, 10000)
chains <- ifelse(short_run, 1, 2)

params_inc <- params_real$inc_par2
si_vec <- seq(min_valid_si, max_valid_si)
## For s3/s4 mix
## cowling_data <- data_s3_s4mix
s3data <- list(
  N = nrow(cowling_data), si = cowling_data$si, max_shed = max_shed,
  alpha2 = params_inc[["shape"]], beta2 = 1 / params_inc[["scale"]],
  M = length(si_vec), si_vec = si_vec, width = width
)
s3pairs <- list(
  N = nrow(data_discrete_pairs), si = data_discrete_pairs$si, max_shed = max_shed,
  alpha2 = params_inc[["shape"]], beta2 = 1 / params_inc[["scale"]],
  M = length(si_vec), si_vec = si_vec, width = width
)
s3s4mix <- list(
  N = nrow(data_s3_s4mix), si = data_s3_s4mix$si, max_shed = max_shed,
  alpha2 = params_inc[["shape"]], beta2 = 1 / params_inc[["scale"]],
  M = length(si_vec), si_vec = si_vec, width = width
)
s3s4mixdiscrete <- list(
  N = nrow(data_discrete_pairs_s3_s4mix), si = data_discrete_pairs_s3_s4mix$si, max_shed = max_shed,
  alpha2 = params_inc[["shape"]], beta2 = 1 / params_inc[["scale"]],
  M = length(si_vec), si_vec = si_vec, width = width
)

#######

fit_model <- function(mixture, recall, right_bias, model_prefix, standata = s3data, obs = cowling_data) {
  prefix <- glue("{model_prefix}_skew_normal")
  infile <- glue("stan-models/{prefix}.stan")
  message(infile)
  if (!file.exists(infile)) message("Does not exist ", infile)
  if(mixture) {
    standata$max_invalid_si <- max_invalid_si
    standata$min_invalid_si <- min_invalid_si
  }
  if (right_bias| recall) standata$nu <- obs$nu
  fit <- stan(
    file = infile, data = standata,  verbose = FALSE, iter = iter,
    chains = chains
  )
  fit
}

fit_model_to_pairs  <- function(mixture, recall, right_bias, model_prefix, standata = s3pairs, obs = data_discrete_pairs) {
  prefix <- glue("{model_prefix}_skew_normal")
  infile <- glue("stan-models/{prefix}.stan")
  message(infile)
  if (!file.exists(infile)) message("Does not exist ", infile)
  if(mixture) {
    standata$max_invalid_si <- max_invalid_si
    standata$min_invalid_si <- min_invalid_si
  }
  if (right_bias| recall) standata$nu <- obs$nu
  fit <- stan(
    file = infile, data = standata,  verbose = FALSE, iter = iter,
    chains = chains
  )
  fit
}


s3s4_model <- function(mixture, recall, right_bias, model_prefix, standata = s3s4mix, obs = data_s3_s4mix) {
  prefix <- glue("{model_prefix}_skew_normal")
  infile <- glue("stan-models/{prefix}.stan")
  message(infile)
  if (!file.exists(infile)) message("Does not exist ", infile)
  if(mixture) {
    standata$max_invalid_si <- max_invalid_si
    standata$min_invalid_si <- min_invalid_si
  }
  if (right_bias| recall) standata$nu <- obs$nu
  fit <- stan(
    file = infile, data = standata,  verbose = FALSE, iter = iter,
    chains = chains
  )
  fit
}


fit_model_to_s3s4pairs  <- function(mixture, recall, right_bias, model_prefix, standata = s3pairs, obs = s3s4mixdiscrete) {
  prefix <- glue("{model_prefix}_skew_normal")
  infile <- glue("stan-models/{prefix}.stan")
  message(infile)
  if (!file.exists(infile)) message("Does not exist ", infile)
  if(mixture) {
    standata$max_invalid_si <- max_invalid_si
    standata$min_invalid_si <- min_invalid_si
  }
  if (right_bias| recall) standata$nu <- obs$nu
  fit <- stan(
    file = infile, data = standata,  verbose = FALSE, iter = iter,
    chains = chains
  )
  fit
}
