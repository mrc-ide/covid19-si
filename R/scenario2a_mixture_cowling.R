data <- readRDS("data/cowling_data_clean.rds")

data_pos <- data %>%
  filter(onset_first_iso > 0) %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso)

width <- min(data_pos$si[data_pos$si > 0]) / 2

fit_mixture_pos <- stan(
  file = here::here("stan-models/scenario2a_mixture_general.stan"),
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
    min_si = min(data_pos$si) - 0.001,
    width = width
  ),
  chains = 2,
  iter = 5000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)

best_params <- map_estimates(fit_mixture)

si_posterior <- simulate_2a_mix(
  mean_inc_og, sd_inc_og,
  list(shape1 = best_params[["alpha1"]],
       shape2 = best_params[["beta1"]]),
  mean_iso, sd_iso,
  best_params[["pinvalid"]], 10000
)

si_posterior <- si_posterior[[3]]

out <- quantile_as_df(si_posterior$si)

p2 <- ggplot() +
  geom_density(
    aes(data_pos$si, fill = "blue"), alpha = 0.3, binwidth = 1
    ) +
  geom_density(
    aes(si_posterior$si, fill = "red"), alpha = 0.3, binwidth = 1
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  theme_minimal() +
  theme(legend.title = element_blank())



ggsave("figures/posterior_serial_interval_2a_mix_cowling_all.png", p2)
