max_shed <- 21
mean_inf <- 4
sd_inf <- 2
mean_inc <- 2
sd_inc <- 1
mean_iso <- 2
sd_iso <- 2
param_offset <- -1
recall <- 5

params_inf <- beta_muvar2shape1shape2(
  (mean_inf - offset) / (max_shed - offset),
  sd_inf^2 /(max_shed - offset)^2
)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = sd_inc/ mean_inc
)

alpha2 <- params_inc$shape
beta2 <- 1 / params_inc$scale

params_iso <- epitrix::gamma_mucv2shapescale(
  mu = mean_iso, cv = sd_iso / mean_iso
)

nsim_post_filter <- 1000

sim_data <- better_simulate_si(
  params_inc, params_inf, params_iso, offset, max_shed, 1e4
)

unfiltered <- sim_data
sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]



sim_data$p_si <- exp(-recall * abs(sim_data$si - sim_data$nu))

idx <- sample(
  nrow(sim_data), nsim_post_filter, replace = TRUE, sim_data$p_si
)
with_recall_bias <- sim_data[idx, ]

ggplot() +
  geom_density(
    aes(unfiltered$t_1, fill = "blue"), col = NA, alpha = 0.2
  ) +
  geom_density(
    aes(sim_data$t_1, fill = "red"), col = NA, alpha = 0.2
  ) +
  geom_density(
    aes(with_recall_bias$t_1, fill = "green"), col = NA, alpha = 0.2
  ) +
  scale_fill_identity(
    breaks = c("blue", "red", "green"),
    labels = c("All", "Filtered", "With recall bias"),
    guide = "legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())


ggplot() +
  geom_point(data = sim_data, aes(nu, si, alpha = p_si)) +
  theme_minimal()

max_si <- ceiling(max(unfiltered$si))
y_vec <- seq(offset, max_si, by = 1)
si_vec <- seq(offset + 0.1 + 0.001, max_si, 1)
width <- 0.1

fit <- stan(
  file = here::here("stan-models/full_model.stan"),
  data = list(
    N = length(with_recall_bias$si),
    si = with_recall_bias$si,
    nu = with_recall_bias$nu,
    max_shed = max_shed,
    offset1 = offset,
    max_si = max(with_recall_bias$si) + 0.001,
    min_si = offset, ## assuming the smallest incubation period is 0
    alpha2 = alpha2,
    beta2 = beta2,
    width = width,
    M = length(si_vec),
    y_vec = si_vec
  ),
  chains = 2,
  iter = 100,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)



## offset should be negative
posterior_inf_params <- function(fit, nsim, max_shed, offset) {
  params <- rstan::extract(fit)
  idx <- sample(length(params[[1]]), nsim, replace = TRUE)
  shape1 <- params[["alpha1"]][idx]
  shape2 <- params[["beta1"]][idx]
  out <- beta_shape1shape22muvar(shape1, shape2)
  out[["mu"]] <- max_shed * out[["mu"]] + offset
  out[["sigma2"]] <- max_shed * max_shed * out[["sigma2"]]
  out
}

out <- posterior_inf_params(fits[[1]], 5000, max_shed, -1)


params_inf_true <- params[[param_grid$params_inf[1]]]

p1 <- ggplot() +
  geom_density(aes(out[[1]]), fill = "red", col = NA, alpha = 0.3) +
  geom_vline(
    xintercept = params_inf_true$mean_inf
  ) +
  geom_vline(
    xintercept = c(
      mean(out[[1]]),
      mean(out[[1]]) - sd(out[[1]]),
      mean(out[[1]]) + sd(out[[1]])
    ), linetype = "dashed"
  ) +
  expand_limits(x = 0) +
  theme_minimal() +
  xlab("Mean infectious period")


p2 <- ggplot() +
  geom_density(aes(sqrt(out[[2]])), fill = "red", col = NA, alpha = 0.3) +
  geom_vline(
    xintercept = params_inf_true$sd_inf
  ) +
  geom_vline(
    xintercept = c(
      mean(sqrt(out[[2]])),
      mean(sqrt(out[[2]])) - sd(sqrt(out[[2]])),
      mean(sqrt(out[[2]])) + sd(sqrt(out[[2]]))
    ), linetype = "dashed"
  ) +
  expand_limits(x = 0) +
  theme_minimal() +
  xlab("SD infectious period")

p <- cowplot::plot_grid(p1, p2, ncol = 1)
