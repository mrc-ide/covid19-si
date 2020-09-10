simulated_1 <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  mean_inf = mean_inf,
  sd_inf = sd_inf,
  nsim = 100
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
  chains = 2,
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
  geom_histogram(alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed")

p2 <- ggplot(NULL, aes(x[["sd"]])) +
  geom_histogram(alpha = 0.3) +
  geom_vline(xintercept = sd_inf, linetype = "dashed")

p <- cowplot::plot_grid(p1, p2, ncol = 1)

ggsave("infectious_profile_params.png", p)


## Simulate with draws from posterior
## nsamples <- length(fitted_params[[1]])
## idx <- sample(nsamples, size = 3000, replace = TRUE)
## ## The two parameters are highly correlated. So we want to sample
## ## them as tuples rather than independently.
## shape <- fitted_params[[1]][idx]
## rate <- fitted_params[[2]][idx]
## posterior_params <- epitrix::gamma_shapescale2mucv(shape, 1 / rate)
## simulated_post <- map2_dfr(
##   posterior_params[[1]],
##   posterior_params[[2]],
##   function(mu, cv) simulate_si(mean_inc, sd_inc, mu,  cv * mu, nsim = 20)
## )

## psi <- ggplot() +
##   geom_histogram(
##     data = simulated_1, aes(x = si, y = ..density..),
##     alpha = 0.3, fill = "blue"
##   ) +

##   geom_histogram(
##     data = simulated_post, aes(x = si, y = ..density..),
##     alpha = 0.3, fill = "red"
##   ) +
##   geom_vline(
##     xintercept = mean(simulated_1$si), col = "red", linetype = "dashed"
##   ) +
##   theme_minimal() +
##   xlab("Serial Interval")

## ggsave("posterior_serial_interval.png", psi)









