## Set up params grid
param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2", "inf_par3"),
  params_inc = c("inc_par1", "inc_par2"),
  params_pinv = c("pinvalid1", "pinvalid2", "pinvalid3"),
  stringsAsFactors = FALSE
)

nsim <- 500
alpha_invalid <- 0.5
beta_invalid <- 0.5
min_si <- -2

plot_1a_mix <- function(valid_si, invalid_si, simulated_si) {
  p1 <- ggplot() +
    geom_histogram(aes(valid_si$si, fill = "blue"), alpha = 0.3) +
    geom_histogram(aes(invalid_si, fill = "red"), alpha = 0.3) +
    theme_minimal() +
    scale_fill_identity(
      guide = "legend",
      breaks = c("blue", "red"),
      labels = c("Valid", "Invalid")
    ) +
    xlim(0, NA) +
    xlab("Serial Interval") +
    theme(legend.title = element_blank())

  p2 <- ggplot() +
    geom_histogram(aes(simulated_si), fill = "purple", alpha = 0.3) +
    theme_minimal() +
    xlim(0, NA) +
    xlab("Serial Interval")


  p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "hv")

  p
}

simulate_1a_mix <- function(mean_inc, sd_inc, params_inf, pinvalid, nsim) {
  valid_si <- simulate_si(
    mean_inc, sd_inc, params_inf$shape1, params_inf$shape2,
    max_shed, NULL, NULL, nsim = nsim
  )

  max_si <- max(valid_si$si)

  invalid_si <- (max_si - min_si)*
    rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

  invalid_si <- invalid_si + min_si
  toss <- runif(nsim)

  simulated_si <- c(
    invalid_si[toss <= pinvalid], valid_si$si[toss > pinvalid]
  )

  list(
    valid_si = valid_si,
    invalid_si = invalid_si,
    simulated_si = simulated_si
  )
}

map_estimates <- function(fit) {

  fitted_params <- rstan::extract(fit)
  map_idx <- which.max(fitted_params[["lp__"]])
  map(fitted_params, ~ .[map_idx])
}

simulated_data <- pmap(
  param_grid,
  function(params_inf, params_inc, params_pinv) {
    outfile <- paste(
      params_inf, params_inc, params_pinv, sep = "_"
    )

    params_inf <- params[[params_inf]]
    params_inc <- params[[params_inc]]
    pinvalid <- params[[params_pinv]]
    params_inf <- beta_muvar2shape1shape2(
      params_inf$mean_inf/max_shed, params_inf$sd_inf^2 / max_shed^2
    )
    simulated_si <- simulate_1a_mix(
      params_inc$mean_inc, params_inc$sd_inc, params_inf, pinvalid, nsim
    )

    p <- plot_1a_mix(
      simulated_si[[1]], simulated_si[[2]], simulated_si[[3]]
    )
    ggsave(glue::glue("figures/{outfile}.png"), p)

    simulated_si <- simulated_si[[3]]
    width <- min(simulated_si[simulated_si > 0]) / 2
    params2_inc <- gamma_mucv2shapescale(
      params_inc$mean_inc, params_inc$sd_inc/ params_inc$mean_inc
    )
    fit_mixture <- stan(
      file = here::here("stan-models/scenario1a_mixture_general.stan"),
      data = list(
        N = length(simulated_si),
        si = simulated_si,
        max_shed = max_shed,
        alpha2 = params2_inc[["shape"]],
        beta2 = 1 / params2_inc[["scale"]],
        alpha_invalid = alpha_invalid,
        beta_invalid = beta_invalid,
        max_si = max(simulated_si) + 0.001,
        min_si = min(simulated_si) - 0.001,
        width = width
      ),
      chains = 5,
      iter = 2000,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
    saveRDS(fit_mixture, glue::glue("stanfits/{outfile}.rds"))

    best_params <- map_estimates(fit_mixture)
    posterior_si <- simulate_1a_mix(
      params_inc$mean_inc, params_inc$sd_inc,
      list(shape1 = best_params[["alpha1"]],
           shape2 = best_params[["beta1"]]),
      best_params[["pinvalid"]], 10000
    )
    saveRDS(
      posterior_si, glue::glue("stanfits/posterior_si_{outfile}.rds")
    )
  }
)

