# run scenario 1 on the cowling data

# load the data

data <- readRDS("data/cowling_data_clean.rds")
data_pos <- data%>%
  filter(si>0)%>%
  filter(onset_first_iso>0)%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)

# restricted
data_pos_restricted <- data_pos%>%
  filter(si<20)

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

# simulate the SI from the mode of the posterior




## p1 <- ggplot(NULL, aes(fitted_params[[2]])) +
##   geom_density(alpha = 0.3, fill = "red") +
##   geom_vline(xintercept = params_inf[[2]], linetype = "dashed") +
##   theme_minimal()

## p2 <- ggplot(NULL, aes(x[["sd"]])) +
##   geom_density(alpha = 0.3, fill = "red") +
##   geom_vline(xintercept = sd_inf, linetype = "dashed") +
##   theme_minimal()

# plot incubation period and empirical SI


x <- seq(0,25, by = 0.01)
y <- dgamma(x = x, shape = params_inc$shape, scale = params_inc$scale)
incub <- data.frame(x,y)

mean_inc <- 5.1 #days
sd_inc <- 3.94 #days
shape_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = (sd_inc/mean_inc)
)[["shape"]]
rate_inc <- 1/epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = (sd_inc/mean_inc)
)[["scale"]]

y_og <- dgamma(x = x, shape = shape_inc, rate = rate_inc)
incub_og <- data.frame(x,y_og)


ggplot(data = data_pos)+
  geom_histogram(aes(x = si, y = ..density..), binwidth = 1)+
  geom_line(data = incub, aes(x = x, y = y))+
  geom_line(data = incub_og, aes(x = x, y = y_og), col = "red")+
  xlab("days")+
  theme_bw()


