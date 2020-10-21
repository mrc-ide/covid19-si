recall_true <- 0.3
nsim <- 500
simulated_2 <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  shape1_inf = params_inf$shape1,
  shape2_inf = params_inf$shape2,
  max_shed = max_shed,
  mean_iso = mean_iso,
  sd_iso = sd_iso,
  nsim = nsim
)

simulated_2 <- simulated_2[simulated_2$t_1 < simulated_2$nu, ]
##simulated_2$p_si <- exp(-recall_true * abs(simulated_2$si - simulated_2$nu))
simulated_2$p_si <- exp(
  abs(simulated_2$t_1 - simulated_2$nu) * -recall_true
)
## Sample with probability p_si
idx <- sample(nrow(simulated_2), nsim, replace = TRUE, prob = simulated_2$p_si)

sampled <- simulated_2[idx, ]

p <- ggplot() +
  geom_point(data = simulated_2[-idx, ], aes(nu, si, col = "black")) +
  geom_point(data = simulated_2[idx, ], aes(nu, si, col = "red")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  expand_limits(x = 0, y = 0) +
  scale_color_identity(
    guide = "legend",
    breaks = c("red", "black"),
    labels = c("Sampled", "Not Sampled")
  ) +
  theme_minimal() +
  xlab("Delay to isolation") +
  ylab("Serial Interval") +
  theme(legend.position = "bottom", legend.title = element_blank())


ggsave("figures/simulated_isolation_si.png", p)


p <- ggplot() +
  geom_histogram(
    data = simulated_2, aes(si, fill = "black"),
    alpha = 0.5, binwidth = 1
  ) +
  geom_histogram(
    data = simulated_2[idx, ], aes(si, fill = "red"),
    alpha = 0.2, binwidth = 1
  ) +
  scale_fill_identity(
    guide = "legend",
    breaks = c("red", "black"),
    labels = c("Sampled", "Simulated")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.position = "bottom", legend.title = element_blank())


ggsave("figures/sampled_with_recall_bias.png", p)


fits_2a_recall <- stan(
  file = here::here("stan-models/scenario2a_recall_bias.stan"),
  data = list(
    N = nrow(sampled),
    si = sampled$si,
    nu = sampled$nu,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    width = min(sampled$si) / 2
  ),
  chains = 3,
  iter = 1000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fitted_params <- rstan::extract(fits_2a_recall)

idx <- which.max(fitted_params[["lp__"]])
shape1 <- fitted_params[["alpha1"]][idx]
shape2 <- fitted_params[["beta1"]][idx]
recall <- fitted_params[["recall"]][idx]


xposterior <- simulate_si(
  mean_inc, sd_inc, shape1, shape2, max_shed, mean_iso, sd_iso, nsim = nsim
)
xposterior <- xposterior[xposterior$t_1 < xposterior$nu, ]

xposterior$p_si <- exp(
  abs(xposterior$t_1 - xposterior$nu) * -recall
)

idx <- sample(
  nrow(xposterior), nrow(xposterior), replace = TRUE,
  prob = xposterior$p_si
)

xposterior <- xposterior [idx, ]

p2 <- ggplot() +
  geom_density(aes(sampled$si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(xposterior$si, fill = "red"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
