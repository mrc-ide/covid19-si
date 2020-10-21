alpha_invalid <- 0.5
beta_invalid <- 0.5

data <- readRDS("data/cowling_data_clean.rds")

data_pos <- data %>%
  filter(si > 0) %>%
  filter(onset_first_iso > 0) %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso)


fit_mixture_pos <- stan(
  file = here::here("stan-models/scenario2a_mixture.stan"),
  data = list(
    N = length(data_pos$si),
    si = data_pos$si,
    nu = data_pos$nu,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data_pos$si) + 0.001,
    min_si = min(data_pos$si),
    width = min(data_pos$si) / 2
  ),
  chains = 2,
  iter = 5000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)


fitted_params <- rstan::extract(fit_mixture_pos)
idx <- which.max(fitted_params[["lp__"]])
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
pinv <- fitted_params[["pinvalid"]][idx]


nsim <- nrow(data_pos)
xposterior <- simulate_si(
  mean_inc, sd_inc, shape1, shape2, max_shed, NULL, NULL, nsim
)


invalid_si <- max(data_pos$si) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

toss <- runif(nsim)

si_posterior <- c(
  invalid_si[toss <= pinv], xposterior[toss > pinv, "si"]
)



p2 <- ggplot() +
  geom_histogram(
    aes(data_pos$si, fill = "blue"), alpha = 0.3, binwidth = 1
    ) +
  geom_histogram(
    aes(si_posterior, fill = "red"), alpha = 0.3, binwidth = 1
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  theme_minimal() +
  theme(legend.title = element_blank())



ggsave("figures/posterior_serial_interval_2a_mix_cowling_positive.png", p2)
