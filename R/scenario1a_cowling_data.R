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
test_fit <- ggmcmc(ggs(fits_1a), here::here("figures/1a_500.pdf"))

## extract fits to turn alpha and beta into mu and cv

fitted_params <- rstan::extract(fits_1a)

max_index <- which(fitted_params$lp__==max(fitted_params$lp__)) 

fitted_max <- c(alpha1 = fitted_params$alpha1[max_index], beta1 = fitted_params$beta1[max_index])

x <- max_shed *
  rbeta(
    n = 100000,
    shape1 = fitted_max[["alpha1"]],
    shape2 = fitted_max[["beta1"]]
  )

p1 <- ggplot() +
  geom_density(aes(x, fill = "forestgreen"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  xlim(0,40)+
  theme_minimal()


ggsave("figures/infectious_profile_params_1a.png", p1)

## simulate si from the most likely incubation period distibution
shape1_max <- fitted_max["alpha1"]
shape2_max <- fitted_max["beta1"]
si_post_max <- simulate_si(mean_inc, sd_inc, shape1_max, shape2_max, max_shed, 2, 2, nsim = 100000)

psi <- ggplot() +
  geom_density(
    data = data_pos, aes(si, fill = "blue"),
    alpha = 0.3
  ) +
  
  geom_density(
    data = si_post_max, aes(si, fill = "red"),
    alpha = 0.3
  ) +
  # geom_density(aes(x, fill = "black"),
  #   alpha = 0.3
  # ) +
  geom_vline(
    xintercept = mean(data_pos$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post_max$si), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())+
  xlim(0,40)+
  ylim(0,0.15)

ggsave("figures/SI_1a.png", psi)

##################################
# we now need to get the 95% CrI #
##################################

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[[1]][idx]
shape2 <- fitted_params[[2]][idx]

# 2. for each infectious profile, simulate an SI distribution
si_post <- matrix(nrow = 1000, ncol = length(idx))
for(i in 1:length(idx)){
si_post[,i] <- (simulate_si(mean_inc, sd_inc, shape1[i], shape2[i], max_shed, 2, 2, nsim = 1000))$si
}

# 3. for each simulated SI distribution, extract the median, mean and sd
si_medians <- colMedians(si_post)
si_means <- colMeans(si_post)
si_sd <- colSds(si_post)

# 4. plot the distribution of each stat and find the 2.5th and 97.5th quantiles
si_mean_95 <- quantile(si_means, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_means), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("mean SI (days)")+
  geom_vline(xintercept = si_mean_95[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_mean_95[2], col = "red", lty = 2)+
  geom_vline(xintercept = mean(si_post_max$si))+
  xlim(3.4, 6.25)
  
si_median_95 <- quantile(si_medians, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(si_post_max$si))+
  xlim(3.4, 6.25)

si_sd_95 <- quantile(si_sd, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(si_post_max$si))+
  xlim(3.4, 6.25)


# 7. plot the simulates SIs as density plots

si_post_df <- data.frame(si_post)
si_post_long <- melt(si_post_df)

p3 <- ggplot(si_post_long, aes(group = variable))+
  geom_density(aes(value), fill = "red", color = NA, alpha = 0.01)+
  theme_bw()+
  geom_density(aes(si), data = si_post_max, colour = "white", inherit.aes = FALSE)


# plot incubation period and empirical SI

# this was to help diagnose why beta wasn't fitting well

x <- seq(0,25, by = 0.01)

shape_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = (sd_inc/mean_inc)
)[["shape"]]
rate_inc <- 1/epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = (sd_inc/mean_inc)
)[["scale"]]

y_og <- dgamma(x = x, shape = shape_inc, rate = rate_inc)
incub_og <- data.frame(x,y_og)


pinc <- ggplot(data = data_pos)+
  geom_histogram(aes(x = si, y = ..density..), binwidth = 1)+
  geom_line(data = incub, aes(x = x, y = y))+
  geom_line(data = incub_og, aes(x = x, y = y_og), col = "red")+
  xlab("days")+
  theme_bw()

ggsave("figures/inc_1a.png", pinc)
