data <- readRDS("data/cowling_data_clean.rds")

data_pos <- data %>%
  filter(si > 0) %>%
  filter(onset_first_iso > 0) %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso)

fits_2a_recall <- stan(
  file = here::here("stan-models/scenario2a_recall_bias.stan"),
  data = list(
    N = nrow(data_pos),
    si = data_pos$si,
    nu = data_pos$nu,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    width = min(data_pos$si) / 2,
    max_si = max(data_pos$si)
  ),
  chains = 3,
  iter = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fitted_params <- rstan::extract(fits_2a_recall)

idx <- which.max(fitted_params[["lp__"]])
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
recall <- fitted_params[["recall"]][idx]

nsim <- 1000

xposterior <- simulate_si(
  mean_inc, sd_inc, shape1, shape2, max_shed, mean_iso, sd_iso, nsim = nsim
)
xposterior <- xposterior[xposterior$t_1 < xposterior$nu, ]

xposterior$p_si <- exp(
  abs(xposterior$si - xposterior$nu) * -recall
)

idx <- sample(
  nrow(xposterior), nrow(xposterior), replace = TRUE,
  prob = xposterior$p_si
)

xposterior_sampled <- xposterior[idx, ]

p2 <- ggplot() +
  geom_density(aes(data_pos$si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(xposterior$si, fill = "black"), alpha = 0.3) +
  geom_density(aes(xposterior_sampled$si, fill = "red"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior", "Posterior Sampled"),
    breaks = c("blue", "black", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave("figures/posterior_si_2a_with_recall_cowling.png", p2)
