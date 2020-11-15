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



sim_data <- readRDS("data/simulated_1.rds")
sim_data <- sim_data[sim_data$type == "valid filtered", ]
sim_data <- sim_data[sample(nrow(sim_data), nsim), ]
width <- min(sim_data$si[sim_data$si > 0]) / 2
fits_4a <- stan(
  file = here::here("stan-models/scenario4a.stan"),
  data = list(
    N = nrow(sim_data),
    si = sim_data$si,
    nu = sim_data$nu,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    offset1 = offset,
    width = 0.1
  ),
  chains = 3,
  iter = 2000,
  verbose = TRUE
)

## out <- pmap(
##   sim_data,
##   function(type, si, nu) {
##     if (si < 0 | nu < 0) return(NULL)
##     else {
##       scenario4a_lpdf(
##         si, nu, max_shed, offset, params_inf[[1]], params_inf[[2]],
##         params_inc[["shape"]], 1 / params_inc[["scale"]], width
##       )
##     }
##   }
## )

best_params <- map_estimates(fits_4a)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = params_inc[[1]], cv = params_inc[[2]] / params_inc[[1]]
)
params_inc$rate <- 1 / params_inc$scale

posterior_si <- better_simulate_si(
  params_inc,
  list(shape1 = best_params[["alpha1"]], shape2 = best_params[["beta1"]]),
  params_iso, offset, max_shed, 10000
)

ggplot() +
  geom_density(
    data = sim_data, aes(si, fill = "red"), alpha = 0.3, col = NA
  ) +
  geom_density(
    data = posterior_si, aes(si, fill = "blue"),
    alpha = 0.3, col = NA
  ) +
  scale_fill_identity(
    breaks = c("red", "blue"),
    labels = c("Simulated", "Posterior"),
    guide  = "legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

true_params <- beta_muvar2shape1shape2(
  params[["inf_par1"]][["mean_inf"]] / max_shed,
  params[["inf_par1"]][["sd_inf"]]^2 / max_shed^2
)

inf_draws <- max_shed *
  rbeta( 10000, true_params$shape1, true_params$shape2)

ggplot() +
  geom_density(
    data = NULL, aes(inf_draws, fill = "red"),
    alpha = 0.3, col = NA) +
  geom_density(
    data = posterior_si, aes(t_1, fill = "blue"), alpha = 0.3, col = NA
  ) +
  scale_fill_identity(
    breaks = c("red", "blue"),
    labels = c("Simulated", "Posterior"),
    guide  = "legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

