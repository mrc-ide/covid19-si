set.seed(42)
nsim <- 1000
alpha_invalid <- 1
beta_invalid <- 1

data <- readRDS("data/cowling_data_clean.rds")

data <- data %>%
  filter(!is.na(onset_first_iso))%>%
  mutate(si = as.numeric(si))%>%
  mutate(nu = as.numeric(onset_first_iso))

data_pos <- data%>%
  filter(si>0)%>%
  filter(onset_first_iso>0)%>%
  mutate(si = as.numeric(si))%>%
  mutate(nu = as.numeric(onset_first_iso))

## fitting the mixture model to positive data
fit_mixture <- stan(
  file = here::here("stan-models/scenario1a_mixture.stan"),
  data = list(
    N = length(data_pos$si),
    si = data_pos$si,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]], #used the original incubation period as it works better with the real data
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data_pos$si) + 0.001,
    width = min(data_pos$si) / 2
  ),
  chains = 5,
  iter = 2000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)

# fitting the mixture model to the whole dataset
fit_mixture <- stan(
  file = here::here("stan-models/scenario1a_mixture_general.stan"),
  data = list(
    N = length(data$si),
    si = data$si,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data$si) + 0.001,
    min_si = min(data$si),
    width = min(data_pos$si) / 2
  ),
  chains = 1,
  iter = 2000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)

test_fit <- ggmcmc(ggs(fit_mixture), here::here("figures/1a_mix_data.pdf"))

fitted_params <- rstan::extract(fit_mixture)
map_idx <- which.max(fitted_params[["lp__"]])
map_params <- map(fitted_params, ~ .[map_idx])
## Simulate with draws from posterior
##nsamples <- length(fitted_params[["alpha1"]])
##idx <- sample(nsamples, size = ceiling(nsamples / 2), replace = FALSE)
shape1 <- map_params[["alpha1"]]
shape2 <- map_params[["beta1"]]
pinv <- map_params[["pinvalid"]]

xposterior <- simulate_si(
  mean_inc_og, sd_inc_og, shape1, shape2, max_shed, nsim = nsim
)
xposterior$si <- round(xposterior$si)

invalid_si <- max(xposterior$si) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

invalid_si <- round(invalid_si)

toss <- runif(nsim)

simulated_si_post <- c(
  invalid_si[toss <= pinv], xposterior$si[toss > pinv]
)


p1 <- ggplot() +
  geom_density(aes(xposterior$t_1, fill = "red"), alpha = 0.3) +
  xlab("Infectious profile") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/posterior_infectious_profile_1a_mix.png", p1)

p2 <- ggplot() +
  geom_density(aes(data_pos$si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(simulated_si_post, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean(data_pos$si), linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())



ggsave("figures/posterior_serial_interval_1a_mix.png", p2)


fitted_params <- rstan::extract(fit_mixture, permuted = FALSE, inc_warmup = TRUE)
pinv_chains <- fitted_params[, , 1]
alpha_chains <- fitted_params[, , 2]
beta_chains <- fitted_params[, , 3]


ggplot() +
  geom_point(aes(alpha_chains[, 1], beta_chains[, 1]), col = "red") +
  geom_point(aes(alpha_chains[, 2], beta_chains[, 2]), col = "blue") +
  geom_point(aes(alpha_chains[, 3], beta_chains[, 3]), col = "green") +
  geom_point(aes(alpha_chains[, 4], beta_chains[, 4]), col = "yellow") +
  geom_point(aes(alpha_chains[, 5], beta_chains[, 5])) +
  theme_minimal() +
  xlab("alpha1") +
  ylab("beta1")


ggplot() +
  geom_point(aes(pinv_chains[, 1], beta_chains[, 1]), col = "red") +
  geom_point(aes(pinv_chains[, 2], beta_chains[, 2]), col = "blue") +
  geom_point(aes(pinv_chains[, 3], beta_chains[, 3]), col = "green") +
  geom_point(aes(pinv_chains[, 4], beta_chains[, 4]), col = "yellow") +
  geom_point(aes(pinv_chains[, 5], beta_chains[, 5])) +
  theme_minimal() +
  xlab("pinv") +
  ylab("beta1")




