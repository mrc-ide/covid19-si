set.seed(42)
nsim <- 2000
alpha_invalid <- 1
beta_invalid <- 1
pinvalid <- 0.1

valid_si <- simulate_si(
  mean_inc, sd_inc, params_inf$shape1, params_inf$shape2,
  max_shed, nsim = nsim
)
valid_si$nu <- rgamma(
  nsim, shape = params_iso$shape, scale = params_iso$scale
)

valid_si <- valid_si[valid_si$t_1 < valid_si$nu, ]

nsim <- nrow(valid_si)

invalid_si <- max(valid_si$si) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

## invalid_inc <- rgamma(
##   nsim, shape = params_inc$shape, scale = params_inc$scale
## )

## invalid_inf <- invalid_si - invalid_inc

invalid_iso <- rgamma(
  nsim, shape = params_iso$shape, scale = params_iso$scale
)

invalid_si = data.frame(si = invalid_si, nu = invalid_iso)

p1 <- ggplot() +
  geom_histogram(aes(valid_si$si, fill = "blue"), alpha = 0.3) +
  geom_histogram(aes(invalid_si$si, fill = "red"), alpha = 0.3) +
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
simulated_si <- rbind(
  invalid_si[toss <= pinvalid, ],
  valid_si[toss > pinvalid, c("si", "nu")]
)



p2 <- ggplot() +
  geom_histogram(aes(simulated_si$si), fill = "purple", alpha = 0.3) +
  theme_minimal() +
  xlim(0, NA) +
  xlab("Serial Interval")

p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "hv")

ggsave("figures/simulated_data_2a_mix.png", p)

fit_mixture <- stan(
  file = here::here("stan-models/scenario2a_mixture.stan"),
  data = list(
    N = nrow(simulated_si),
    si = simulated_si$si,
    nu = simulated_si$nu,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    width = min(simulated_si$si) / 2,
    max_si = max(simulated_si$si) + 0.001
  ),
  chains = 5,
  iter = 5000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)


fitted_params <- rstan::extract(fit_mixture)

## Simulate with draws from posterior
##nsamples <- length(fitted_params[["alpha1"]])
##idx <- sample(nsamples, size = ceiling(nsamples / 2), replace = FALSE)
idx <- which.max(fitted_params[["lp__"]])
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
pinv <- fitted_params[["pinvalid"]][idx]


xposterior <- simulate_si(
  mean_inc, sd_inc, shape1, shape2, max_shed, nsim = nsim
)

toss <- runif(nsim)

si_posterior <- rbind(
  invalid_si[toss <= pinv, ],
  xposterior[toss > pinv, "si"]
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


ggsave("figures/posterior_infectious_profile_2a_mix.png", p1)

p2 <- ggplot() +
  geom_density(aes(simulated_si$si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(si_posterior$si, fill = "red"), alpha = 0.3) +
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



ggsave("figures/posterior_serial_interval_2a_mix.png", p2)


fitted_params <- rstan::extract(fit_mixture, permuted = FALSE, inc_warmup = TRUE)
pinv_chains <- fitted_params[, , "pinvalid"]
alpha_chains <- fitted_params[, , "alpha1"]
beta_chains <- fitted_params[, , "beta1"]


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

