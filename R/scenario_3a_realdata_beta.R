## s3 + mixture REAL DATA + Beta distribution for
## for inf profile
####################
# read in the data #
####################

data <- readRDS("data/cowling_data_clean.rds")

real_data <- data %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso) %>%
  dplyr::filter(!is.na(nu))

real_data <- real_data[order(real_data$nu),]



###################################################
# select assumed parameters (from list in global) #
###################################################

params_offset <- params_real$offset3
params_inc <- params_real$inc_par2
max_shed <- params_real$maxshed3
first_valid_nu <- match(params_offset + 1, real_data$nu)
######################
# fit the stan model #
######################
si_vec <- seq(params_offset + 1, max_si)


fit_3a_real <- stan(
  file = here::here("stan-models/scenario3a_mixture_beta.stan"),
  data = list(
    N = length(real_data$si),
    si = real_data$si,
    max_shed = max_shed,
    offset1 = params_offset,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_invalid_si = max_si,
    min_invalid_si = min_invalid_si,
    width = width,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
  ),
  ##chains = 1, iter = 100,
  seed = 42,
  verbose = TRUE
)

saveRDS(fit_3a_real, "stanfits/release/scenario3a_mixture_beta.rds")
