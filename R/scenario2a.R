simulated_2 <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  shape1_inf = params_inf$shape1,
  shape2_inf = params_inf$shape2,
  max_shed = max_shed,
  mean_iso = 2,
  sd_iso = 2,
  nsim = 1000
)

simulated_2 <- simulated_2[simulated_2$t_1 < simulated_2$nu, ]

## grid <- expand.grid(
##   alpha1 = seq(1, 100, by = 1), beta1 = seq(0.1, 100, by = 0.5)
## )

## out <- pmap(
##   grid,
##   function(alpha1, beta1) {
##     y <- g(simulated_2$si, simulated_2$nu, alpha1, beta1,
##       params_inc[["shape"]], 1 / params_inc[["scale"]]
##       )
##     sum(y)
##   }
## )



fits_2a <- stan(
  file = here::here("stan-models/scenario2a.stan"),
  data = list(
    N = nrow(simulated_2),
    si = simulated_2$si,
    nu = simulated_2$nu,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]]
  ),
  chains = 3,
  iter = 10000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fitted_params <- rstan::extract(fits_2a)

x <- hermione::beta_shape1shape22muvar(
  fitted_params[["alpha1"]], fitted_params[["beta1"]]
)

x[["mu"]] <- max_shed * x[["mu"]]
x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
x[["sd"]] <- sqrt(x[["sigma2"]])

p1 <- ggplot(NULL, aes(x[["mu"]])) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed")

p2 <- ggplot(NULL, aes(x[["sd"]])) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = sd_inf, linetype = "dashed")

p <- cowplot::plot_grid(p1, p2, ncol = 1)

ggsave("infectious_profile_params.png", p)

nsamples <- length(fitted_params[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[[1]][idx]
shape2 <- fitted_params[[2]][idx]

si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2)

psi <- ggplot() +
  geom_density(
    data = simulated_2, aes(si),
    alpha = 0.3, fill = "blue"
  ) +

  geom_density(
    data = si_post, aes(si),
    alpha = 0.3, fill = "red"
  ) +
  geom_vline(
    xintercept = mean(simulated_2$si), col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Serial Interval")

ggsave("posterior_serial_interval.png", psi)
