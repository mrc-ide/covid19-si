# s3 + mixture REAL DATA

# read in the data

data <- readRDS("data/cowling_data_clean.rds")

real_data <- data%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)%>%
  dplyr::filter(!is.na(nu))

# select assumed parameters (from list in global)

params_offset <- params_real$offset3
params_inc <- params_real$inc_par1
max_shed <- params_real$maxshed3

# fit the stan model 

fit_3a_real_3 <- stan(
  file = here::here("stan-models/scenario3a_mixture_general.stan"),
  data = list(
    N = length(real_data$si),
    si = real_data$si,
    max_shed = max_shed,
    offset1 = params_offset,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(real_data$si) + 0.001,
    min_si = min(real_data$si) - 0.001,
    width = width
  ),
  chains = 2,
  iter = 2000,
  seed = 42,
  verbose = TRUE
)

 v # check mcmc diagnostics
diagnos <- ggmcmc(ggs(fit_3a_real_3), here::here("3a_real_3.pdf"))

# 