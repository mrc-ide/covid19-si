valid_si <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  shape1_inf = params_inf$shape1,
  shape2_inf = params_inf$shape2,
  max_shed = max_shed,
  nsim = 500
)
invalid_si <- runif(500, min(valid_si$si), max(valid_si$si))
pinvalid <- 0.3

toss <- runif(500)

simulated_si <- c(
  invalid_si[toss <= pinvalid], valid_si$si[toss > pinvalid]
)

ggplot() + geom_density(aes(simulated_si))
## simulated_si <- simulated_si[1:10]


fit_mixture  <- stan(
  file = here::here("stan-models/scenario1a_mixture.stan"),
  data = list(
    N = length(simulated_si),
    si = simulated_si,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = 1,
    beta_invalid = 1,
    max_si = max(simulated_si),
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
pinv <- fitted_params[["pinvalid"]][idx]

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
  theme(legend.title = element_blank())

ggsave("posterior_distr_pinvalid.png", p)

params <- data.frame(shape1 = shape1, shape2 = shape2, p = pinv)
## g <- Vectorize(simulate_mixture)

si_post <- pmap_dfr(
  params,
  function(shape1, shape2, p) {
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
