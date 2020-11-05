# run scenario 3 + mixture model on the cowling data

# offset
offset <- 3

# select number of simulations from posterior, and parameters for invalid dist.
nsim <- 10000
alpha_invalid <- 0.5
beta_invalid <- 0.5

# load the data
data <- readRDS("data/cowling_data_clean.rds")
data <- data%>%
    mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)

# sub-set to only incude those SIs that are possible under our assumed offset
data_offset <- data%>%
  filter(si>-offset)%>%
  filter(nu>-offset)%>%


# fit the model
fits_3a_mix <- stan(
  file = here::here("stan-models/scenario3a_mixture_general.stan"),
  data = list(
    N = nrow(data),
    si = data$si,
    max_shed = 21,
    offset = offset,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data$si) + 0.001,
    min_si = min(data$si) - 0.001
  ),
  chains = 1,
  iter = 3000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

test_fit_3a_mix <- ggmcmc(ggs(fits_3a_mix), here::here("figures/3a_mix.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params_3a_mix <- rstan::extract(fits_3a_mix)
saveRDS(fitted_params_3a_mix, file = "fitted_params_3a_mix.rds")
max_index <- which(fitted_params_3a_mix$lp__==max(fitted_params_3a_mix$lp__)) 

fitted_max <- c(alpha1 = fitted_params_3a_mix$alpha1[max_index], beta1 = fitted_params_3a_mix$beta1[max_index])

## x <- hermione::beta_shape1shape22muvar(
##   fitted_params_3a_mix[["alpha1"]], fitted_params_3a_mix[["beta1"]]
## )

## x[["mu"]] <- max_shed * x[["mu"]]
## x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
## x[["sd"]] <- sqrt(x[["sigma2"]])

x <- (max_shed *
        rbeta(
          n = 100000,
          shape1 = fitted_max[["alpha1"]],
          shape2 = fitted_max[["beta1"]]
        ))-offset

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

incubation <- rgamma(
  n = 100000,
  shape = params_inc_og[["shape"]],
  scale = params_inc_og[["scale"]]
)

si_post <- x + incubation

psi <- ggplot() +
  geom_histogram(
    data = data_offset, aes(si, y = ..density.., fill = "blue"),
    alpha = 0.3,
    binwidth = 1
  ) +
  
  geom_density(aes(si_post, fill = "red"),
               alpha = 0.3
  ) +
  # geom_density(aes(x, fill = "black"),
  #   alpha = 0.3
  # ) +
  geom_vline(
    xintercept = mean(data_offset$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())
ggsave("figures/SI_3a_mix.png", psi)

##################################
# we now need to get the 95% CrI #
##################################

si_post_max3a_mix <- si_post

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params_3a_mix[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params_3a_mix[[1]][idx]
shape2 <- fitted_params_3a_mix[[2]][idx]

# 2. for each infectious profile, simulate an SI distribution (of 1000 SIs)
si_post_3a_mix <- matrix(nrow = 1000, ncol = length(idx))
for(i in 1:length(idx)){
  si_post_3a_mix[,i] <- (simulate_si(mean_inc_og, sd_inc_og, shape1[i], shape2[i], max_shed, 2, 2, nsim = 1000))$si - offset
}

# 3. for each simulated SI distribution, extract the median, mean and sd
library(matrixStats)
si_medians_3a_mix <- colMedians(si_post_3a_mix)
si_means_3a_mix <- colMeans(si_post_3a_mix)
si_sd_3a_mix <- colSds(si_post_3a_mix)

# 4. plot the distribution of each stat and find the 2.5th and 97.5th quantiles
si_mean_95_3a_mix <- quantile(si_means_3a_mix, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_means_3a_mix), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("mean SI (days)")+
  geom_vline(xintercept = si_mean_95_3a_mix[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_mean_95_3a_mix[2], col = "red", lty = 2)+
  geom_vline(xintercept = mean(si_post_max3a_mix))+
  xlim(4, 20)

si_median_95_3a_mix <- quantile(si_medians_3a_mix, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians_3a_mix), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95_3a_mix[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95_3a_mix[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(si_post_max3a_mix))+
  xlim(4, 20)

si_sd_95_3a_mix <- quantile(si_sd_3a_mix, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd_3a_mix), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95_3a_mix[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95_3a_mix[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(si_post_max3a_mix))+
  xlim(4, 20)