# run scenario 3 + mixture model on the cowling data

# offset
offset <- -3

# select number of simulations from posterior, and parameters for invalid dist.
nsim <- 10000
alpha_invalid <- 1
beta_invalid <- 1

# load the data
data <- readRDS("data/cowling_data_clean.rds")

data_c <- data%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)%>%
  dplyr::filter(!is.na(nu))

# sub-set to only incude those SIs that are possible under our assumed offset
data_offset <- data_c %>%
  filter(si > -offset) %>%
  filter(nu > -offset)


# fit the model
fits_4a_mix <- stan(
  file = here::here("stan-models/scenario4a_mixture.stan"),
  data = list(
    N = nrow(data_c),
    si = data_c$si,
    nu = data_c$nu,
    max_shed = 21,
    offset1 = offset,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data_c$si) + 0.001,
    min_si = min(data_c$si) - 0.001,
    width = 0.1
  ),
  chains = 4,
  iter = 5000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

test_fit_4a_mix <- ggmcmc(ggs(fits_4a_mix), here::here("figures/4a_mix_uniform.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params_4a_mix <- rstan::extract(fits_4a_mix)
saveRDS(fitted_params_4a_mix, file = "fitted_params_4a_mix.rds")
max_index <- which(fitted_params_4a_mix$lp__==max(fitted_params_4a_mix$lp__))

fitted_max <- c(alpha1 = fitted_params_4a_mix$alpha1[max_index], 
                beta1 = fitted_params_4a_mix$beta1[max_index], 
                p_invalid = fitted_params_4a_mix$pinvalid[max_index])
shape1_max <- fitted_max["alpha1"]
shape2_max <- fitted_max["beta1"]

x <- (max_shed *
        rbeta(
          n = 100000,
          shape1 = shape1_max,
          shape2 = shape2_max
        ))+offset

p1 <- ggplot() +
  geom_density(aes(x, fill = "blue"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  theme_minimal()


ggsave("figures/infectious_profile_params_4a_mix.png", p1)

## simulate si from the most likely incubation period distibution


si_post <- (simulate_si(mean_inc_og, sd_inc_og, shape1_max, shape2_max, max_shed, nsim = 100000))$si + offset

psi <- ggplot() +
  geom_histogram(
    data = data_c, aes(si, y = ..density.., fill = "blue"),
    alpha = 0.3,
    binwidth = 1
  ) +
  
  geom_density(aes(si_post, fill = "red"),
               alpha = 0.3
  ) +
  geom_vline(
    xintercept = mean(data_c$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior (valid SIs)"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  xlim(NA,40)+
  ylim(0,0.105)+
  theme(legend.title = element_blank())
ggsave("figures/SI_4a_mix.png", psi, width = 7, height = 7, units = "in", dpi = 300, device = "png")

## including invalid SIs in the figure

si_post_p <- (simulate_3a_mix(mean_inc_og, sd_inc_og, shape1_max, shape2_max, max_shed,
                              pinvalid = p_invalid_max, nsim = 100000, offset = offset,
                              alpha_invalid, beta_invalid, min_si = min(data_c$si), max_si = max(data_c$si)))
si_post_iv <- si_post_p$simulated_si$si

psi_iv <- ggplot() +
  geom_histogram(
    data = data_c, aes(si, y = ..density.., fill = "blue"),
    alpha = 0.3,
    binwidth = 1
  ) +
  
  geom_density(aes(si_post_iv, fill = "red"),
               alpha = 0.3
  ) +
  geom_vline(
    xintercept = mean(data_c$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post_iv), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior (valid & invalid SIs)"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())
ggsave("figures/SI_4a_mix_invalid.png", psi_iv)


##################################
# we now need to get the 95% CrI #
##################################

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params_4a_mix[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params_4a_mix$alpha1[idx]
shape2 <- fitted_params_4a_mix$beta1[idx]

# 2. for each infectious profile, simulate an SI distribution (of 10000 SIs)
si_post_4a <- matrix(nrow = 1000, ncol = length(idx))
for(i in 1:length(idx)){
  si_post_4a[,i] <- (simulate_si(mean_inc_og, sd_inc_og, shape1[i], shape2[i], max_shed, nsim = 1000))$si + offset
}

# 3. for each simulated SI distribution, extract the median, mean and sd
library(matrixStats)
si_medians_4a <- colMedians(si_post_4a)
si_means_4a <- colMeans(si_post_4a)
si_sd_4a <- colSds(si_post_4a)

# 4. plot the distribution of each stat and find the 2.5th and 97.5th quantiles
si_mean_95_4a <- quantile(si_means_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_means_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("mean SI (days)")+
  geom_vline(xintercept = si_mean_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_mean_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = mean(si_post))

si_median_95_4a <- quantile(si_medians_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(si_post))

si_sd_95_4a <- quantile(si_sd_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(si_post))

## getting the CrI for p_invalid

p_invalid_max <- fitted_max["p_invalid"]

p_invalid_post <- fitted_params_4a_mix$pinvalid

mean(p_invalid_post)
median(p_invalid_post)
quantile(p_invalid_post, c(0.025, 0.975))

## adding conditional fitted distribution to the plot


fitted_params_4a_mix <- readRDS("fitted_params_4a_mix.rds")

inf_dist <- (max_shed *
               rbeta(
                 n = 10000,
                 shape1 = fitted_max[["alpha1"]],
                 shape2 = fitted_max[["beta1"]]
               ))+offset

n <- length(data_c$nu)
inf_filt <- list() 

for(v in 1:n){ #for each observed nu, filter out any sampled inf delays < nu
  
  inf_filt[[v]] <- inf_dist[which(inf_dist<data_c$nu[v])]
}

inf_conditional <- unlist(inf_filt) 
inc <- rgamma(
  n = length(inf_conditional),
  shape = params_inc_og[["shape"]],
  scale = params_inc_og[["scale"]]
)
si_conditional <- inf_conditional + inc

psi+
  geom_density(aes(si_conditional, fill = "green"),
               alpha = 0.3
  )+
  geom_vline(
    xintercept = mean(si_conditional), col = "green", linetype = "dashed"
  )+
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior", "Conditional"),
    breaks = c("blue", "red", "green")
  )

## and with the invalid sis as well to assess fit of the whole model
min_si <- min(data_c$si)
max_si <- max(data_c$si)
nsim <- length(si_conditional)
pinvalid <- p_invalid_max

invalid_si <- (max_si - min_si)*
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

invalid_si <- invalid_si + min_si
toss <- runif(nsim)

sim_invalid <- data.frame(si = invalid_si[toss <= pinvalid], group = "invalid")
sim_valid <- data.frame(si = si_conditional[toss > pinvalid], group = "valid")

conditional_si <- rbind(sim_invalid, sim_valid) 

psi_cond <- psi+
  geom_density(aes(conditional_si$si, fill = "green"),
               alpha = 0.3
  )+
  geom_vline(
    xintercept = mean(conditional_si$si), col = "green", linetype = "dashed"
  )+
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior", "Conditional"),
    breaks = c("blue", "red", "green")
  )
