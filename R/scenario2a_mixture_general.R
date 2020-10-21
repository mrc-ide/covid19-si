alpha_invalid <- 0.5
beta_invalid <- 0.5

data <- readRDS("data/cowling_data_clean.rds")

data_all <- data %>%
  filter(onset_first_iso > 0) %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso)

width <- min(data_all$si[data_all$si > 0]) / 2

fit_mixture <- stan(
  file = here::here("stan-models/scenario2a_mixture_general.stan"),
  data = list(
    N = length(data_all$si),
    si = data_all$si,
    nu = data_all$nu,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data_all$si) + 0.001,
    min_si = min(data_all$si) - 0.001,
    width = width
  ),
  chains = 2,
  iter = 2000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)

fitted_params <- rstan::extract(fit_mixture)
idx <- which.max(fitted_params[["lp__"]])
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
pinv <- fitted_params[["pinvalid"]][idx]


nsim <- nrow(data_all)
xposterior <- simulate_si(
  mean_inc, sd_inc, shape1, shape2, max_shed, NULL, NULL, nsim
)


invalid_si <- (max(data_all$si) - min(data_all$si)) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

invalid_si <- invalid_si + min(data_all$si)
toss <- runif(nsim)

si_posterior <- c(
  invalid_si[toss <= pinv], xposterior[toss > pinv, "si"]
)



p2 <- ggplot() +
  geom_histogram(aes(data_all$si, fill = "blue"), alpha = 0.3) +
  geom_histogram(aes(si_posterior, fill = "red"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())



ggsave("figures/posterior_serial_interval_2a_mix_cowling_all.png", p2)
