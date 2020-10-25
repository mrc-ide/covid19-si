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
