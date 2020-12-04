# run scenario 3 on the cowling data

# offset
offset <- -3

# load the data
data <- readRDS("data/cowling_data_clean.rds")

# sub-set to only incude those SIs that are possible under our assumed offset
data_offset <- data%>%
  filter(si>offset)%>%
  filter(onset_first_iso>offset)%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)

# fit the model
fits_4a <- stan(
  file = here::here("stan-models/scenario4a.stan"),
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    nu = data_offset$nu,
    max_shed = 21,
    offset1 = offset,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
  width = 0.1),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


## Check convergence etc, using ggmcmc
test_fit <- ggmcmc(ggs(fits_4a), here::here("figures/4a.pdf"))

## extract fits to turn alpha and beta into mu and cv


fitted_params <- rstan::extract(fits_4a)
saveRDS(fitted_params, file = "fitted_params_4a.rds")
#fitted_params <- readRDS("fitted_params_4a.rds")
max_index <- which(fitted_params$lp__==max(fitted_params$lp__)) 

fitted_max <- c(alpha1 = fitted_params$alpha1[max_index], beta1 = fitted_params$beta1[max_index])

## x <- hermione::beta_shape1shape22muvar(
##   fitted_params[["alpha1"]], fitted_params[["beta1"]]
## )

## x[["mu"]] <- max_shed * x[["mu"]]
## x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
## x[["sd"]] <- sqrt(x[["sigma2"]])

x <- ((max_shed-offset) *
        rbeta(
          n = 100000,
          shape1 = fitted_max[["alpha1"]],
          shape2 = fitted_max[["beta1"]]
        ))+offset

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

ggsave("figures/SI_4a.png", psi)

##################################
# we now need to get the 95% CrI #
##################################

fitted_params_4a <- fitted_params
si_post_max4a <- si_post

# 1. sample alpha and beta from the posterior to get the infectious profile

nsamples <- length(fitted_params_4a[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params_4a[[1]][idx]
shape2 <- fitted_params_4a[[2]][idx]

# 2. for each infectious profile, simulate an SI distribution (of 1000 SIs)
si_post_4a <- matrix(nrow = 10000, ncol = length(idx))
for(i in 1:length(idx)){
  params_inf <- list(shape1 = shape1[i], shape2 = shape2[i])
  si_post_4a[,i] <- better_simulate_si(params_inc = params_inc_og, params_inf = params_inf, 
                                       params_iso = params_iso, min_inf = offset, max_inf = max_shed, 
                                       nsim = 10000)$si}

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
  geom_vline(xintercept = mean(si_post_max4a))

si_median_95_4a <- quantile(si_medians_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_medians_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("median SI (days)")+
  geom_vline(xintercept = si_median_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_median_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = median(si_post_max4a))

si_sd_95_4a <- quantile(si_sd_4a, c(0.025, 0.975))
ggplot()+
  geom_density(aes(si_sd_4a), fill = "red", alpha = 0.3)+
  theme_bw()+
  xlab("sd SI (days)")+
  geom_vline(xintercept = si_sd_95_4a[1], col = "red", lty = 2)+
  geom_vline(xintercept = si_sd_95_4a[2], col = "red", lty = 2)+
  geom_vline(xintercept = sd(si_post_max4a))

## adding conditional fitted distribution to the plot

  inf_dist <- (max_shed *
                 rbeta(
                   n = 10000,
                   shape1 = fitted_max[["alpha1"]],
                   shape2 = fitted_max[["beta1"]]
                 ))+offset
  
  n <- length(data_offset$nu)
  inf_filt <- list() 
  
  for(v in 1:n){ #for each observed nu, filter out any sampled inf delays < nu

  inf_filt[[v]] <- inf_dist[which(inf_dist<data_offset$nu[v])]
  }
 
  inf_conditional <- unlist(inf_filt) 
  inc <- rgamma(
    n = length(inf_conditional),
    shape = params_inc_og[["shape"]],
    scale = params_inc_og[["scale"]]
  )
  si_conditional <- inf_conditional + inc

  psi_cond <- psi+
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
  ggsave("figures/SI_4a.png", psi_cond, width = 7, height = 7, units = "in", dpi = 300, device = "png")
  
  