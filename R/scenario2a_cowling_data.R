# load the data

data <- readRDS("data/cowling_data_clean.rds")
data_pos <- data%>%
  filter(si>0)%>%
  filter(onset_first_iso>0)%>%
  mutate(si = as.numeric(si))%>%
  mutate(nu = as.numeric(onset_first_iso))

data_pos_test <- data_pos%>%
  filter(nu<21)


fits_2a <- stan(
  file = here::here("stan-models/scenario2a.stan"),
  data = list(
    N = nrow(data_pos),
    si = data_pos$si,
    nu = data_pos$nu,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]]
  ),
  chains = 1,
  iter = 1000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fitted_params <- rstan::extract(fits_2a)

x <- beta_shape1shape22muvar(
  fitted_params[["alpha1"]], fitted_params[["beta1"]]
)

x[["mu"]] <- max_shed * x[["mu"]]
x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
x[["sd"]] <- sqrt(x[["sigma2"]])

p1 <- ggplot(NULL, aes(x[["mu"]])) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed")

p2 <- ggplot(NULL, aes(x[["sd"]])) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = sd_inf, linetype = "dashed")

p <- cowplot::plot_grid(p1, p2, ncol = 1)

ggsave("infectious_profile_params.png", p)

nsamples <- length(fitted_params[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[[1]][idx]
shape2 <- fitted_params[[2]][idx]

si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2)

psi <- ggplot() +
  geom_density(
    data = simulated_2, aes(si),
    alpha = 0.3, fill = "blue"
  ) +
  
  geom_density(
    data = si_post, aes(si),
    alpha = 0.3, fill = "red"
  ) +
  geom_vline(
    xintercept = mean(simulated_2$si), col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Serial Interval")

ggsave("posterior_serial_interval.png", psi)
