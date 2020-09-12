simulated_1 <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  shape1_inf = params_inf$shape1,
  shape2_inf = params_inf$shape2,
  max_shed = max_shed,
  nsim = 500
)

fits_1a <- stan(
  file = here::here("stan-models/scenario1a.stan"),
  data = list(
    N = nrow(simulated_1),
    si = simulated_1$si,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]]
  ),
  chains = 3,
  iter = 5000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


## Check convergence etc, using ggmcmc
## test_fit <- ggmcmc(ggs(fit_1a), here::here("figures/1a.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params <- rstan::extract(fits_1a)

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

ggsave("figures/infectious_profile_params_1a.png", p)


## Simulate with draws from posterior
nsamples <- length(fitted_params[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[[1]][idx]
shape2 <- fitted_params[[2]][idx]

si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2)

psi <- ggplot() +
  geom_density(
    data = simulated_1, aes(si),
    alpha = 0.3, fill = "blue"
  ) +

  geom_density(
    data = si_post, aes(si),
    alpha = 0.3, fill = "red"
  ) +
  geom_vline(
    xintercept = mean(simulated_1$si), col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Serial Interval")

ggsave("figures/posterior_serial_interval_1a.png", psi)






## p1 <- ggplot(NULL, aes(fitted_params[[2]])) +
##   geom_density(alpha = 0.3, fill = "red") +
##   geom_vline(xintercept = params_inf[[2]], linetype = "dashed") +
##   theme_minimal()

## p2 <- ggplot(NULL, aes(x[["sd"]])) +
##   geom_density(alpha = 0.3, fill = "red") +
##   geom_vline(xintercept = sd_inf, linetype = "dashed") +
##   theme_minimal()




