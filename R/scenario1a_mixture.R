simulate_mixture <- function(params_ip, params_inf, params_invalid, pinvalid,
                             max_shed = 21, nsim = 500) {


  t_1 <- max_shed *
    rbeta(
      n = nsim, shape1 = params_inf$shape1, shape2 = params_inf$shape2
    )

  t_2 <- stats::rgamma(
    n = nsim, shape = params_ip$shape, scale = params_ip$scale
    )

  valid_si <- t_1 + t_2
  invalid_t1 <- max_shed * rbeta(
    nsim, params_invalid$shape1, params_invalid$shape2
  )
  invalid_si <- invalid_t1 + t_2
  toss <- runif(nsim)
  c(invalid_si[toss <= pinvalid], valid_si[toss > pinvalid])
}

pinvalid <- 0.3
params_invalid <- list(shape1 = 0.8, shape2 = 0.8)

simulated_si <- simulate_mixture(
  params_inc, params_inf, params_invalid, pinvalid, nsim = 100
)
## Need to make sure SI is less than max_shed otherwise stan cribs
simulated_si <- simulated_si[simulated_si < 21]


##simulated_si <- simulated_si[1:10]


fit_mixture <- stan(
  file = here::here("stan-models/scenario1a_mixture.stan"),
  data = list(
    N = length(simulated_si),
    si = simulated_si,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = params_invalid$shape1,
    beta_invalid = params_invalid$shape2
  ),
  chains = 3,
  iter = 5000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


fitted_params <- rstan::extract(fit_mixture)

x <- hermione::beta_shape1shape22muvar(
  fitted_params[["alpha1"]], fitted_params[["beta1"]]
)

x[["mu"]] <- max_shed * x[["mu"]]
x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
x[["sd"]] <- sqrt(x[["sigma2"]])

p1 <- ggplot(NULL, aes(x[["mu"]])) +
  geom_density(alpha = 0.3, fill = "red") +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  theme_minimal()

p2 <- ggplot(NULL, aes(x[["sd"]])) +
  geom_density(alpha = 0.3, fill = "red") +
  geom_vline(xintercept = sd_inf, linetype = "dashed") +
  theme_minimal()

p <- cowplot::plot_grid(p1, p2, ncol = 1)

ggsave("figures/infectious_profile_params_1a_mix.png", p)


## Simulate with draws from posterior
nsamples <- length(fitted_params[["alpha1"]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
pinv <- fitted_params[["theta"]][idx, 1]

params <- data.frame(shape1 = shape1, shape2 = shape2, p = pinv)
##g <- Vectorize(simulate_mixture)

si_post <- pmap_dfr(
  params,
  function(shape1, shape2, p){
    out <- simulate_mixture(
      params_inc,
      params_inf = list(shape1 = shape1, shape2 = shape2),
      params_invalid = params_invalid,
      pinvalid = p,
      max_shed = max_shed,
      nsim = 10
    )
    data.frame(si = out)
  }
)

psi <- ggplot() +
  geom_density(
    data = NULL, aes(simulated_si),
    alpha = 0.3, fill = "blue"
  ) +

  geom_density(
    data = si_post, aes(si),
    alpha = 0.3, fill = "red"
  ) +
  geom_vline(
    xintercept = mean(simulated_si), col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Serial Interval")

ggsave("figures/posterior_serial_interval_1a_mix.png", psi)


