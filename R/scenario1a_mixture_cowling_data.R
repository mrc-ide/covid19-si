set.seed(42)
nsim <- 1000
alpha_invalid <- 0.5
beta_invalid <- 0.5

# run scenario 1 mixture model on the cowling data

# load the data

simulated_si <- readRDS("data/cowling_data_clean.rds") %>%
  filter(onset_first_iso > 0) %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso)


width <- min(simulated_si$si[simulated_si$si > 0]) / 2

fit_mixture <- stan(
  file = here::here("stan-models/scenario1a_mixture_general.stan"),
  data = list(
    N = length(simulated_si$si),
    si = simulated_si$si,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(simulated_si$si) + 0.001,
    min_si = min(simulated_si$si) - 0.001,
    width = width
  ),
  chains = 2,
  iter = 5000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)
min_si <- min(simulated_si$si)
best_params <- map_estimates(fit_mixture)
posterior_si <- simulate_1a_mix(
  mean_inc_og, sd_inc_og,
  list(shape1 = best_params[["alpha1"]],
       shape2 = best_params[["beta1"]]),
  best_params[["pinvalid"]], 10000
)


posterior_si <- posterior_si[[3]]

p2 <- ggplot() +
  geom_density(aes(simulated_si$si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(posterior_si, fill = "red"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

out <- quantile_as_df(posterior_si)


ggsave("figures/posterior_serial_interval_1a_mix_cowling.png", p2)
