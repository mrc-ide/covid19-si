sim_data <- readRDS("data/simulated_1.rds")
sim_data <- sim_data[sim_data$type == "mixed and recall", ]

## nsim is the number of data points to be fed to stan
nsim <- 250

sim_data <- sim_data[sample(nrow(sim_data), nsim), ]

params_inc <- gamma_mucv2shapescale(
  mu = params[["inc_par1"]][["mean_inc"]],
  cv = params[["inc_par1"]][["sd_inc"]] /
    params[["inc_par1"]][["mean_inc"]]
)

width <- min(sim_data$si[sim_data$si > 0]) / 2

fits_1a <- stan(
  file = here::here("stan-models/scenario1a.stan"),
  data = list(
    N = nrow(sim_data),
    si = sim_data$si,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    width = width
  ),
  chains = 3,
  iter = 5000,
  verbose = TRUE
)

best_params <- map_estimates(fits_1a)

posterior_si <- simulate_si(
  params[["inc_par1"]][["mean_inc"]],
  params[["inc_par1"]][["sd_inc"]],
  shape1 = best_params[["alpha1"]],
  shape2 = best_params[["beta1"]],
  max_shed, NULL, NULL, nsim = 10000
)

ggplot() +
  geom_density(data = sim_data, aes(si), fill = "red", alpha = 0.3) +
  geom_density(data = posterior_si, aes(si), fill = "blue", alpha = 0.3) +
  theme_minimal()

true_params <- beta_muvar2shape1shape2(
  params[["inf_par1"]][["mean_inf"]] / max_shed,
  params[["inf_par1"]][["sd_inf"]]^2 / max_shed^2
)

inf_draws <- max_shed *
  rbeta( 10000, true_params$shape1, true_params$shape2)

ggplot() +
  geom_density(data = NULL, aes(inf_draws), fill = "red", alpha = 0.3) +
  geom_density(data = posterior_si, aes(t_1), fill = "blue", alpha = 0.3) +
  theme_minimal()
