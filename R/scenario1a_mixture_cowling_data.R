# set.seed(42)
# nsim <- 1000
alpha_invalid <- 1
beta_invalid <- 1
# pinvalid <- 0.3
# 
# valid_si <- simulate_si(
#   mean_inc, sd_inc, params_inf$shape1, params_inf$shape2, max_shed, nsim = nsim
# )
# 
# invalid_si <- max(valid_si$si) *
#   rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)
# 
# p1 <- ggplot() +
#   geom_density(aes(valid_si$si, fill = "blue"), alpha = 0.3) +
#   geom_density(aes(invalid_si, fill = "red"), alpha = 0.3) +
#   theme_minimal() +
#   scale_fill_identity(
#     guide = "legend",
#     breaks = c("blue", "red"),
#     labels = c("Valid", "Invalid")
#   ) +
#   xlim(0, NA) +
#   xlab("Serial Interval") +
#   theme(legend.title = element_blank())
# 
# 
# 
# 
# toss <- runif(nsim)
# 
# simulated_si <- c(
#   invalid_si[toss <= pinvalid], valid_si$si[toss > pinvalid]
# )
# 
# p2 <- ggplot() +
#   geom_density(aes(simulated_si), fill = "purple", alpha = 0.3) +
#   theme_minimal() +
#   xlim(0, NA) +
#   xlab("Serial Interval")
# 
# p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "hv")
# 
# ggsave("figures/simulated_data_1a_mix.png", p)

# run scenario 1 mixture model on the cowling data

# load the data

data <- readRDS("data/cowling_data_clean.rds")
data_pos <- data%>%
  filter(si>0)%>%
  filter(onset_first_iso>0)%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)


fit_mixture <- stan(
  file = here::here("stan-models/scenario1a_mixture.stan"),
  data = list(
    N = length(data_pos$si),
    si = data_pos$si,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data_pos$si) + 0.001,
    min_si = min(data_pos$si),
    width = min(data_pos$si) / 2
  ),
  chains = 2,
  iter = 5000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)

test_fit <- ggmcmc(ggs(fit_mixture), here::here("figures/1a_mixture.pdf"))

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
  mean_inc, sd_inc, shape1, shape2, max_shed, nsim = nsim
)


p1 <- ggplot() +
  geom_density(aes(valid_si$t_1, fill = "blue"), alpha = 0.3) +
  geom_density(aes(xposterior$t_1, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Infectious profile") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


ggsave("figures/posterior_infectious_profile_1a_mix.png", p1)

p2 <- ggplot() +
  geom_density(aes(simulated_si, fill = "blue"), alpha = 0.3) +
  geom_density(aes(xposterior$si, fill = "red"), alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed") +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Posterior"),
    breaks = c("blue", "red")
  ) +
  xlab("Serial Interval") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(legend.title = element_blank())



ggsave("figures/posterior_serial_interval_1a_mix.png", p2)
