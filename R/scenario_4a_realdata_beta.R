# s4 + mixture REAL DATA
###################################################
# select assumed parameters (from list in global) #
###################################################

#params_offset <- params_real$offset3
#params_offset <- -6
params_offset <- -11
params_inc <- params_real$inc_par2
max_shed <- params_real$maxshed3

####################
# read in the data #
####################

data <- readRDS("data/cowling_data_clean.rds")

real_data <- data%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)%>%
  dplyr::filter(!is.na(nu))

real_data <- real_data[order(real_data$nu),]

first_valid_nu <- match(params_offset + 2, real_data$nu) #changed to 2 for params_offset = 11 bc no nu = 10


######################
# fit the stan model #
######################
si_vec <- seq(params_offset + 1, max_si)

fit_4a_real_beta_11 <- stan(
  file = here::here("stan-models/scenario4a_mixture_beta.stan"),
  data = list(
    N = length(real_data$si),
    si = real_data$si,
    nu = real_data$nu,
    max_shed = max_shed,
    offset1 = params_offset,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    max_invalid_si = max_si,
    min_invalid_si = min_invalid_si,
    width = width,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu =  first_valid_nu
  ),
  seed = 42,
  verbose = TRUE
)

##########################
# check mcmc diagnostics #
##########################

diagnos <- ggmcmc(ggs(fit_4a_real_beta), here::here("4a_real_beta11.pdf"))
saveRDS(fit_4a_real_beta_11, file = "fit_4a_11.rds")
fit_4a_real_beta <- readRDS("fit_4a_03032021.rds")

################
# extract fits #
################
fitted_params_4a_real_beta <- rstan::extract(fit_4a_real_beta_11)

  # best fits
max_index <- which.max(fitted_params_4a_real_beta$lp__)
params_inf_max <- list(
  shape1 = fitted_params_4a_real_beta$alpha1[max_index],
  shape2 = fitted_params_4a_real_beta$beta1[max_index]
)

p_invalid_max <- fitted_params_4a_real_beta$pinvalid[max_index]

##################
# plot best fits #
##################

  # plot the best fitting infectious profile
shifted_inf <- ((max_shed-params_offset) *
        rbeta(
          n = 100000,
          shape1 = params_inf_max$shape1,
          shape2 = params_inf_max$shape2
        )) + params_offset

p_inf <- ggplot() +
  geom_density(aes(shifted_inf, fill = "green"), alpha = 0.3) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  xlab("delay from symptoms to transmission (days)")+
  theme_minimal()+
  theme(axis.text=element_text(size=15))

p_presymp3 <- sum(shifted_inf>0) / length(shifted_inf)

## Mean and CV of positive delays
mean_delay <- mean(real_data$nu[real_data$nu > params_offset])
sd_delay <- sd(real_data$nu[real_data$nu > params_offset])
iso_params_obs <- epitrix::gamma_mucv2shapescale(mean_delay,
                                                 sd_delay / mean_delay
                                                 )

### plot the resulting SI of best fitting parameters
posterior_si <- better_simulate_si(
  params_inc, params_inf_max, iso_params_obs, params_offset, max_shed,
  nsim_pre_filter
)
posterior_si <- posterior_si[posterior_si$t_1 <= posterior_si$nu, ]

### plot with invalid serial intervals
posterior_tot_si <- simulate_4a_mix(params_inc, params_inf_max, iso_params_obs,
                                    params_offset, max_shed, p_invalid_max, nsim_pre_filter,
                                                alpha_invalid, beta_invalid, min_invalid_si, max_si)
posterior_all_si <- posterior_tot_si$simulated_si
posterior_uncon_si <- posterior_tot_si$unconditional

psi <- ggplot() +
  geom_histogram(
    data = real_data, aes(si, y = ..density.., fill = "blue"),
    alpha = 0.3,
    binwidth = 1
  ) +
  geom_density(
    aes(posterior_si$si, fill = "green"), color = NA, alpha = 0.3
  ) +
  geom_density(
    aes(posterior_uncon_si$si, fill = "darkgreen"), color = NA, alpha = 0.3
  ) +
  geom_density(
    aes(posterior_all_si$si), color = "red", alpha = 0.3
  ) +
  geom_vline(
    xintercept = mean(real_data$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(posterior_si$si), col = "green", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(posterior_uncon_si$si), col = "darkgreen", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior - 'conditional SI'", "Posterior (valid & invalid SIs)"),
    breaks = c("blue", "green", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=15))

###########
# 95% CrI #
###########

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params_4a_real_beta[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params_4a_real_beta$alpha1[idx]
shape2 <- fitted_params_4a_real_beta$beta1[idx]

# 2. for each infectious profile, simulate an SI distribution (of 10000 SIs)
si_post_4a <- matrix(nrow = 10000, ncol = length(idx))
for(i in 1:length(idx)){
  params_inf <- list(shape1 = shape1[i], shape2 = shape2[i])
  si_post_4a[,i] <- (simulate_4a_mix(params_inc = params_inc, params_inf = params_inf,
                                     params_iso = list(shape = 1, scale = 1), offset = params_offset, max_shed,
                                     pinvalid = p_invalid_max, nsim = 10000, alpha_invalid,
                                     beta_invalid, min_si = min(real_data$si), max_si = max(real_data$si)))$valid_si$si
}

# 3. for each simulated SI distribution, extract the median, mean and sd
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
  geom_vline(xintercept = mean(posterior_uncon_si$si))

si_median_95_4a <- quantile(si_medians_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(posterior_uncon_si$si))

si_sd_95_4a <- quantile(si_sd_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(posterior_uncon_si$si))

## getting the CrI for p_invalid

p_invalid_post <- fitted_params_4a_real_beta$pinvalid

mean(p_invalid_post)
median(p_invalid_post)
quantile(p_invalid_post, c(0.025, 0.975))



