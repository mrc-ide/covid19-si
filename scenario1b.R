simulated_1b <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  shape1_inf = params_inf$shape1,
  shape2_inf = params_inf$shape2,
  max_shed = max_shed,
  nsim = 500
)


fits_1b <- stan(
  file = here::here("stan-models/scenario1b.stan"),
  data = list(
    N = nrow(simulated_1b),
    si = simulated_1b$si,
    max_shed = 21,
    width = min(simulated_1b$si) / 2
  ),
  chains = 3,
  iter = 5000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


## Check convergence etc, using ggmcmc
## ggmcmc(ggs(fits_1b), here::here("figures/1b.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params <- rstan::extract(fits_1b)
map_idx <- which.max(fitted_params[["lp__"]])
map_params <- map(fitted_params, ~ .[map_idx])
## x <- hermione::beta_shape1shape22muvar(
##   fitted_params[["alpha1"]], fitted_params[["beta1"]]
## )

## x[["mu"]] <- max_shed * x[["mu"]]
## x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
## x[["sd"]] <- sqrt(x[["sigma2"]])

x1 <- max_shed *
  rbeta(
    n = 10000,
    shape1 = map_params[["alpha1"]],
    shape2 = map_params[["beta1"]]
  )

p1 <- ggplot() +
  geom_density(aes(x1, fill = "blue"), alpha = 0.3) +
  geom_density(aes(simulated_1b$t_1, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Infectious Period") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/infectious_profile_params_1b.png", p1)

x2 <- rgamma(
    n = 10000, shape = map_params[["alpha2"]], rate = map_params[["beta2"]]
  )

p2 <- ggplot() +
  geom_density(aes(x2, fill = "blue"), alpha = 0.3) +
  geom_density(aes(simulated_1b$t_2, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Incubation Period") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/incubation_period_params_1b.png", p2)


x3 <- x1 + x2

p3 <- ggplot() +
  geom_density(aes(x3, fill = "blue"), alpha = 0.3) +
  geom_density(aes(simulated_1b$si, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/serial_interval_params_1b.png", p3)

