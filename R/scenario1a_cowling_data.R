# run scenario 1 on the cowling data

# load the data

data <- readRDS("data/cowling_data_clean.rds")
data_pos <- data%>%
  filter(si>0)%>%
  filter(onset_first_iso>0)%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)

# fit the model

fits_1a <- stan(
  file = here::here("stan-models/scenario1a.stan"),
  data = list(
    N = nrow(data_pos),
    si = data_pos$si,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]]
  ),
  chains = 3,
  iter = 5000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


## Check convergence etc, using ggmcmc
test_fit <- ggmcmc(ggs(fits_1a), here::here("figures/1a.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params <- rstan::extract(fits_1a)

## x <- hermione::beta_shape1shape22muvar(
##   fitted_params[["alpha1"]], fitted_params[["beta1"]]
## )

## x[["mu"]] <- max_shed * x[["mu"]]
## x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
## x[["sd"]] <- sqrt(x[["sigma2"]])

x <- max_shed *
  rbeta(
    n = length(fitted_params[["alpha1"]]),
    shape1 = fitted_params[["alpha1"]],
    shape2 = fitted_params[["beta1"]]
  )

p1 <- ggplot() +
  geom_density(aes(x, fill = "blue"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  theme_minimal()


ggsave("figures/infectious_profile_params_1a.png", p1)


## Simulate with draws from posterior
nsamples <- length(fitted_params[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[[1]][idx]
shape2 <- fitted_params[[2]][idx]

si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2)

psi <- ggplot() +
  geom_density(
    data = data_pos, aes(si, fill = "blue"),
    alpha = 0.3
  ) +
  
  geom_density(
    data = si_post, aes(si, fill = "red"),
    alpha = 0.3
  ) +
  geom_vline(
    xintercept = mean(data_pos$si), col = "blue", linetype = "dashed"
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

ggsave("figures/posterior_serial_interval_1a.png", psi)

# also look at the MLE! 






## p1 <- ggplot(NULL, aes(fitted_params[[2]])) +
##   geom_density(alpha = 0.3, fill = "red") +
##   geom_vline(xintercept = params_inf[[2]], linetype = "dashed") +
##   theme_minimal()

## p2 <- ggplot(NULL, aes(x[["sd"]])) +
##   geom_density(alpha = 0.3, fill = "red") +
##   geom_vline(xintercept = sd_inf, linetype = "dashed") +
##   theme_minimal()




