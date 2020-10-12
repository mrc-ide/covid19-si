# load the data

data <- readRDS("data/cowling_data_clean.rds")
data_pos <- data%>%
  filter(si>0)%>%
  filter(onset_first_iso>0)%>%
  mutate(si = as.numeric(si))%>%
  mutate(nu = as.numeric(onset_first_iso))


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
  chains = 3,
  iter = 5000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)
test_fit <- ggmcmc(ggs(fits_2a), here::here("figures/2a_working.pdf"))
fitted_params_2a <- rstan::extract(fits_2a)

max_index <- which(fitted_params_2a$lp__==max(fitted_params_2a$lp__)) 

fitted_max <- c(alpha1 = fitted_params_2a$alpha1[max_index], beta1 = fitted_params_2a$beta1[max_index])


x <- max_shed *
  rbeta(
    n = 100000,
    shape1 = fitted_max[["alpha1"]],
    shape2 = fitted_max[["beta1"]]
  )

p2 <- ggplot() +
  geom_density(aes(x, fill = "blue"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  theme_minimal()


ggsave("figures/infectious_profile_params_2a.png", p2)

## simulate si from the most likely incubation period distibution
shape1 <- fitted_max["alpha1"]
shape2 <- fitted_max["beta1"]
si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2, nsim = 100000)

p2si <- ggplot() +
  geom_density(
    data = data_pos_test, aes(si, fill = "blue"),
    alpha = 0.3
  ) +
  
  geom_density(
    data = si_post, aes(si, fill = "red"),
    alpha = 0.3
  ) +
  # geom_density(aes(x, fill = "black"),
  #   alpha = 0.3
  # ) +
  geom_vline(
    xintercept = mean(data_pos_test$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post$si), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())

ggsave("figures/SI_2a.png", p2si)


# using the whole posteriors to get the 95% CrI

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
