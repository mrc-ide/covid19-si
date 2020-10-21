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

max_index <- which.max(fitted_params_2a$lp_)

fitted_max <- c(alpha1 = fitted_params_2a$alpha1[max_index],
                beta1 = fitted_params_2a$beta1[max_index])


x <- max_shed *
  rbeta(
    n = 100000,
    shape1 = fitted_max[["alpha1"]],
    shape2 = fitted_max[["beta1"]]
  )

p2 <- ggplot() +
  geom_density(aes(x, fill = "forestgreen"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  xlim(0,40)+
  theme_minimal()


ggsave("figures/infectious_profile_params_2a.png", p2)

## simulate si from the most likely incubation period distibution
shape1 <- fitted_max["alpha1"]
shape2 <- fitted_max["beta1"]
si_post_max2a <- simulate_si(mean_inc, sd_inc, shape1, shape2,
                             max_shed, 2, 2, nsim = 100000)

p2si <- ggplot() +
  geom_density(
    data = data_pos, aes(si, fill = "blue"),
    alpha = 0.3
  ) +

  geom_density(
    data = si_post_max2a, aes(si, fill = "red"),
    alpha = 0.3
  ) +
  # geom_density(aes(x, fill = "black"),
  #   alpha = 0.3
  # ) +
  geom_vline(
    xintercept = mean(data_pos$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post_max2a$si), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  xlim(0,40)+
  ylim(0,0.15)+
  theme(legend.title = element_blank())

ggsave("figures/SI_2a.png", p2si)


# using the whole posteriors to get the 95% CrI

##################################
# we now need to get the 95% CrI #
##################################

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params_2a[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params_2a[[1]][idx]
shape2 <- fitted_params_2a[[2]][idx]

# 2. for each infectious profile, simulate an SI distribution
si_post_2a <- matrix(nrow = 1000, ncol = length(idx))
for(i in 1:length(idx)){
  si_post_2a[,i] <- (simulate_si(mean_inc, sd_inc, shape1[i], shape2[i], max_shed, 2, 2, nsim = 1000))$si
}

# 3. for each simulated SI distribution, extract the median, mean and sd
si_medians_2a <- colMedians(si_post_2a)
si_means_2a <- colMeans(si_post_2a)
si_sd_2a <- colSds(si_post_2a)

# 4. plot the distribution of each stat and find the 2.5th and 97.5th quantiles
si_mean_95_2a <- quantile(si_means_2a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_means_2a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("mean SI (days)")+
  geom_vline(xintercept = si_mean_95_2a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_mean_95_2a[2], col = "red", lty = 2)+
  geom_vline(xintercept = mean(si_post_max2a$si))+
  xlim(4, 20)

si_median_95_2a <- quantile(si_medians_2a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians_2a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95_2a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95_2a[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(si_post_max2a$si))+
  xlim(4, 20)

si_sd_95_2a <- quantile(si_sd_2a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd_2a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95_2a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95_2a[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(si_post_max2a$si))+
  xlim(4, 20)


# 7. plot the simulates SIs as density plots

si_post_2a_df <- data.frame(si_post_2a)
si_post_2a_long <- melt(si_post_2a_df)

p4 <- ggplot(si_post_2a_long, aes(group = variable))+
  geom_density(aes(value), fill = "red", color = NA, alpha = 0.01)+
  theme_bw()+
  geom_density(aes(si), data = si_post_max2a, colour = "white", inherit.aes = FALSE)

