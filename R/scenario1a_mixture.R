
##invalid_si <- runif(100, max(valid_si$si), 3 * max(valid_si$si))

invalid_si <- max(valid_si$si) * rbeta(100, shape1 = 1.5, shape2 = 6.5)

p1 <- ggplot() +
  geom_density(aes(valid_si$si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(invalid_si, fill = "red"), alpha = 0.3) +
  theme_minimal() +
  scale_fill_identity(
    guide = "legend",
    breaks = c("blue", "red"),
    labels = c("Valid", "Invalid")
  ) +
  xlim(0, NA) +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())


pinvalid <- 0.3

toss <- runif(100)

simulated_si <- c(
  invalid_si[toss <= pinvalid], valid_si$si[toss > pinvalid]
)

p2 <- ggplot() +
  geom_density(aes(simulated_si), fill = "purple", alpha = 0.3) +
  theme_minimal() +
  xlim(0, NA) +
  xlab("Serial Interval")

p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "hv")

ggsave("figures/simulated_data_1a_mix.png", p)

fit_mixture  <- stan(
  file = here::here("stan-models/scenario1a_mixture.stan"),
  data = list(
    N = length(simulated_si),
    si = simulated_si,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = 1.5,
    beta_invalid = 6.5,
    max_si = max(simulated_si) + 0.001,
    min_si = min(simulated_si)
  ),
  chains = 2,
  iter = 1000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)


fitted_params <- rstan::extract(fit_mixture)

## Simulate with draws from posterior
nsamples <- length(fitted_params[["alpha1"]])
idx <- sample(nsamples, size = ceiling(nsamples / 2), replace = FALSE)
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
pinv <- sample(fitted_params[["pinvalid"]], ceiling(nsamples / 2))

xvalid <- max_shed *
  rbeta(
    n = idx,
    shape1 = fitted_params[["alpha1"]],
    shape2 = fitted_params[["beta1"]]
  )
xinvalid <- 32 *
  rbeta(
    n = idx,
    shape1 = params_invalid$shape1,
    shape2 = params_invalid$shape2
  )
toss <- runif(idx)
x <- c(
  xinvalid[toss < pinv], xvalid[toss > pinv]
)

p <- ggplot() +
  geom_density(aes(simulated_si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(x, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave("figures/infectious_profile_samples_1a_mix.png", p)



p <- ggplot() +
  geom_density(aes(pinv, fill = "red"), alpha = 0.3) +
  scale_fill_identity(
    guide = "none",
    labels = "Posterior",
    breaks = "red"
  ) +
  xlab("Probability(SI invalid)") +
  geom_vline(xintercept = pinvalid, linetype = "dashed") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  xlim(0, 1)

ggsave("figures/posterior_distr_pinvalid.png", p)

params <- data.frame(shape1 = shape1, shape2 = shape2, p = pinv)
## g <- Vectorize(simulate_mixture)

si_post <- pmap_dfr(
  params,
  function(shape1, shape2, p) {
    valid_si <- simulate_si(
      mean_ip = mean_inc,
      sd_ip = sd_inc,
      shape1_inf = shape1,
      shape2_inf = shape2,
      max_shed = max_shed,
      nsim = 100
    )

    invalid_si <- max(valid_si$si) *
      rbeta(100, shape1 = 1.5, shape2 = 6.5)
    toss <- runif(100)
    out <- c(
      invalid_si[toss <= p], valid_si$si[toss > p]
    )
    data.frame(si = out)
  }
)

psi <- ggplot() +
  geom_density(
    data = NULL, aes(simulated_si, fill = "blue"),
    alpha = 0.3
  ) +
  geom_density(
    data = si_post, aes(si, fill = "red"),
    alpha = 0.3
  ) +
  geom_vline(
    xintercept = mean(simulated_si), col = "red", linetype = "dashed"
  ) +
  xlab("Serial Interval") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/posterior_serial_interval_1a_mix.png", psi)
