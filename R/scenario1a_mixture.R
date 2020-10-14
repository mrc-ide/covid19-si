set.seed(42)
nsim <- 100
alpha_invalid <- 1.5
beta_invalid <- 6.5
pinvalid <- 0.03

valid_si <- simulate_si(
  mean_inc, sd_inc, 6.5, 1.5, max_shed, nsim = nsim
)

##valid_si$si <- round(valid_si$si)

invalid_si <- max(valid_si$si) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

##invalid_si <- round(invalid_si)

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




toss <- runif(nsim)

simulated_si <- c(
  invalid_si[toss <= pinvalid], valid_si$si[toss > pinvalid]
)


p2 <- ggplot() +
  geom_histogram(aes(simulated_si), fill = "purple", alpha = 0.3) +
 ## geom_histogram(aes(simulated_si_post), alpha = 0.3) +
  theme_minimal() +
  xlim(0, NA) +
  xlab("Serial Interval")


p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "hv")

ggsave("figures/simulated_data_1a_mix.png", p)

simulated_si <- simulated_si[simulated_si > 0]
##simulated_si <- simulated_si[!duplicated(simulated_si)]

fit_mixture <- stan(
  file = here::here("stan-models/scenario1a_mixture.stan"),
  data = list(
    N = length(simulated_si),
    si = simulated_si,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(simulated_si) + 0.001,
    width = min(simulated_si) / 2
  ),
  chains = 5,
  iter = 2000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)


fitted_params <- rstan::extract(fit_mixture)
map_idx <- which.max(fitted_params[["lp__"]])
map_params <- map(fitted_params, ~ .[map_idx])
## Simulate with draws from posterior
##nsamples <- length(fitted_params[["alpha1"]])
##idx <- sample(nsamples, size = ceiling(nsamples / 2), replace = FALSE)
shape1 <- map_params[["alpha1"]]
shape2 <- map_params[["beta1"]]
pinv <- map_params[["pinvalid"]]

xposterior <- simulate_si(
  mean_inc, sd_inc, shape1, shape2, max_shed, nsim = nsim
)
xposterior$si <- round(xposterior$si)

invalid_si <- max(xposterior$si) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

invalid_si <- round(invalid_si)

toss <- runif(nsim)

simulated_si_post <- c(
  invalid_si[toss <= pinv], xposterior$si[toss > pinv]
)


p1 <- ggplot() +
  geom_density(aes(valid_si$t_1, fill = "blue"), alpha = 0.3) +
  geom_density(aes(xposterior$t_1, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Infectious profile") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/posterior_infectious_profile_1a_mix.png", p1)

p2 <- ggplot() +
  geom_density(aes(simulated_si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(simulated_si_post, fill = "red"), alpha = 0.3) +
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



ggsave("figures/posterior_serial_interval_1a_mix.png", p2)


fitted_params <- rstan::extract(fit_mixture, permuted = FALSE, inc_warmup = TRUE)
pinv_chains <- fitted_params[, , 1]
alpha_chains <- fitted_params[, , 2]
beta_chains <- fitted_params[, , 3]


ggplot() +
  geom_point(aes(alpha_chains[, 1], beta_chains[, 1]), col = "red") +
  geom_point(aes(alpha_chains[, 2], beta_chains[, 2]), col = "blue") +
  geom_point(aes(alpha_chains[, 3], beta_chains[, 3]), col = "green") +
  geom_point(aes(alpha_chains[, 4], beta_chains[, 4]), col = "yellow") +
  geom_point(aes(alpha_chains[, 5], beta_chains[, 5])) +
  theme_minimal() +
  xlab("alpha1") +
  ylab("beta1")


ggplot() +
  geom_point(aes(pinv_chains[, 1], beta_chains[, 1]), col = "red") +
  geom_point(aes(pinv_chains[, 2], beta_chains[, 2]), col = "blue") +
  geom_point(aes(pinv_chains[, 3], beta_chains[, 3]), col = "green") +
  geom_point(aes(pinv_chains[, 4], beta_chains[, 4]), col = "yellow") +
  geom_point(aes(pinv_chains[, 5], beta_chains[, 5])) +
  theme_minimal() +
  xlab("pinv") +
  ylab("beta1")




