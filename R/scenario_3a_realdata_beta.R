## s3 + mixture REAL DATA + Beta distribution for
## for inf profile
####################
# read in the data #
####################

data <- readRDS("data/cowling_data_clean.rds")

real_data <- data %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso) %>%
  dplyr::filter(!is.na(nu))

real_data <- real_data[order(real_data$nu),]



###################################################
# select assumed parameters (from list in global) #
###################################################

params_offset <- params_real$offset3
params_inc <- params_real$inc_par2
max_shed <- params_real$maxshed3
first_valid_nu <- match(params_offset + 1, real_data$nu)
######################
# fit the stan model #
######################
si_vec <- seq(params_offset + 1, max_si)


fit_3a_real <- stan(
  file = here::here("stan-models/scenario3a_mixture_beta.stan"),
  data = list(
    N = length(real_data$si),
    si = real_data$si,
    max_shed = max_shed,
    offset1 = params_offset,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_invalid_si = max_si,
    min_invalid_si = min_invalid_si,
    width = width,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
  ),
  ##chains = 1, iter = 100,
  seed = 42,
  verbose = TRUE
)


##########################
# check mcmc diagnostics #
##########################

diagnos <- ggmcmc(
  ggs(fit_3a_real), here::here("figures/3a_real.pdf")
)

saveRDS(
  fit_3a_real,
  file = glue::glue("stanfits/fit_3a_{Sys.Date()}.rds")
)

################
# extract fits #
################

fitted_params_3a_real <- rstan::extract(fit_3a_real)
### best fits
max_index <- which.max(fitted_params_3a_real$lp__)
params_inf_max <- list(
  shape1 = fitted_params_3a_real$alpha1[max_index],
  shape2 = fitted_params_3a_real$beta1[max_index]
)
p_invalid_max <- fitted_params_3a_real$pinvalid[max_index]

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
  geom_density(
    aes(shifted_inf, fill = "green"), alpha = 0.3,
    col = NA
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Posterior"),
    breaks = c("red")
  ) +
  xlab("delay from symptoms to transmission (days)")+
  theme_minimal()

p_presymp3 <- sum(shifted_inf > 0) / length(shifted_inf)
  # plot the resulting SI of best fitting parameters
si_post <- (simulate_3a_mix(params_inc = params_inc, params_inf = params_inf_max,
                            params_iso = list(shape = 1, scale = 1), offset = params_offset, max_shed,
                            pinvalid = p_invalid_max, nsim = 5000000, alpha_invalid,
                            beta_invalid, min_si = min(real_data$si), max_si = max(real_data$si)))
# note - in the above, params_iso will not be used and can be any value

si_post_iv <- si_post$simulated_si$si
si_post_v <- si_post$valid_si$si
psi <- ggplot() +
  geom_histogram(
    data = real_data, aes(si, y = ..density.., fill = "blue"),
    alpha = 0.3,
    binwidth = 1
  ) +

  geom_density(aes(si_post_iv, fill = "red"), color = NA,
               alpha = 0.3
  ) +
  geom_density(aes(si_post_v, fill = "green"), color = NA,
               alpha = 0.3
  ) +
  geom_vline(
    xintercept = mean(real_data$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post_iv), col = "red", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_post_v), col = "green", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior - 'true SI'", "Posterior (valid & invalid SIs)"),
    breaks = c("blue", "green", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )


###########
# 95% CrI #
###########

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params_3a_real[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params_3a_real$alpha1[idx]
shape2 <- fitted_params_3a_real$beta1[idx]

# 2. for each infectious profile, simulate an SI distribution (of 10000 SIs)
si_post_3a <- matrix(nrow = 10000, ncol = length(idx))
for(i in 1:length(idx)){
  params_inf <- list(shape1 = shape1[i], shape2 = shape2[i])
  si_post_3a[,i] <- (simulate_3a_mix(params_inc = params_inc, params_inf = params_inf,
                                     params_iso = list(shape = 1, scale = 1), offset = params_offset, max_shed,
                                     pinvalid = p_invalid_max, nsim = 10000, alpha_invalid,
                                     beta_invalid, min_si = min(real_data$si), max_si = max(real_data$si)))$valid_si$si
}

# 3. for each simulated SI distribution, extract the median, mean and sd
si_medians_3a <- colMedians(si_post_3a)
si_means_3a <- colMeans(si_post_3a)
si_sd_3a <- colSds(si_post_3a)

# 4. plot the distribution of each stat and find the 2.5th and 97.5th quantiles
si_mean_95_3a <- quantile(si_means_3a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_means_3a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("mean SI (days)")+
  geom_vline(xintercept = si_mean_95_3a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_mean_95_3a[2], col = "red", lty = 2)+
  geom_vline(xintercept = mean(si_post_v))

si_median_95_3a <- quantile(si_medians_3a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians_3a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95_3a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95_3a[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(si_post_v))

si_sd_95_3a <- quantile(si_sd_3a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd_3a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95_3a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95_3a[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(si_post_v))

## getting the CrI for p_invalid

p_invalid_post <- fitted_params_3a_real$pinvalid

mean(p_invalid_post)
median(p_invalid_post)
quantile(p_invalid_post, c(0.025, 0.975))



