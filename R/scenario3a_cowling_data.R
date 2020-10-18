# run scenario 3 on the cowling data

# set offset (the number of days prior to symptoms that transmission begins)
offset <- 3

# load the data
data <- readRDS("data/cowling_data_clean.rds")

# sub-set to only incude those SIs that are possible under our assumed offset
data_offset <- data%>%
  filter(si>-offset)%>%
  filter(onset_first_iso>-offset)%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)

data_test <- data_offset%>%
  filter(si<0)
# fit the model
fits_3a <- stan(
  file = here::here("stan-models/scenario3a.stan"),
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    max_shed = 21,
    offset = offset,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]]
  ),
  chains = 1,
  iter = 10,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


## Check convergence etc, using ggmcmc
test_fit <- ggmcmc(ggs(fits_1a), here::here("figures/1aogt.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params <- rstan::extract(fits_1a100_rest)

max_index <- which(fitted_params$lp__==max(fitted_params$lp__)) 

fitted_max <- c(alpha1 = fitted_params$alpha1[max_index], beta1 = fitted_params$beta1[max_index])

## x <- hermione::beta_shape1shape22muvar(
##   fitted_params[["alpha1"]], fitted_params[["beta1"]]
## )

## x[["mu"]] <- max_shed * x[["mu"]]
## x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
## x[["sd"]] <- sqrt(x[["sigma2"]])

x <- max_shed *
  rbeta(
    n = 100000,
    shape1 = fitted_max[["alpha1"]],
    shape2 = fitted_max[["beta1"]]
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

## simulate si from the most likely incubation period distibution
shape1 <- fitted_max["alpha1"]
shape2 <- fitted_max["beta1"]
si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2, nsim = 100000)

psi <- ggplot() +
  geom_density(
    data = data_pos, aes(si, fill = "blue"),
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