sim_data <- readRDS("data/simulated_1.rds")
sim_data <- sim_data[sim_data$type == "mixed and recall", ]



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

