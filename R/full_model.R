max_shed <- 21
mean_inf <- 4.973684 # these give an alpha of 4.5 and beta of 14.5: both on the grid for offset = 0
sd_inf <- 1.99637
mean_inf <-4.217391 
mean_sd <- 2.022031
mean_inc <- 3
sd_inc <- 1
## very short isolation
mean_iso <- 2
sd_iso <- 1
offset <- -1

params_inf <- beta_muvar2shape1shape2(
  (mean_inf-offset)/(max_shed-offset), sd_inf^2 /(max_shed-offset)^2
)

params_inf <- list(shape1 = 5, shape2 = 16)

params_inf_mu <- ((max_shed-offset)*beta_shape1shape22muvar(4.5,14.5)$mu)+offset
params_inf_sd <- (max_shed-offset)*sqrt(beta_shape1shape22muvar(4.5,14.5)$sigma2)

params_inf_mu <- ((max_shed-offset)*beta_shape1shape22muvar(5,16)$mu)+offset
params_inf_sd <- (max_shed-offset)*sqrt(beta_shape1shape22muvar(5,16)$sigma2)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = sd_inc/ mean_inc
)

alpha2 <- params_inc$shape
beta2 <- 1 / params_inc$scale

params_iso <- epitrix::gamma_mucv2shapescale(
  mu = mean_iso, cv = sd_iso / mean_iso
)

nsim_post_filter <- 500

expose_stan_functions("stan-models/likelihoods.stan")


mle <- list()
mean <- list()
p <- list()
t_1_posterior<- list()
p2 <- list()
mle_norm <- list()
mean_norm <- list()
p_norm <- list()
t_1_posterior_norm<- list()
p2_norm <- list()
max_si <- vector(length = 10)

sim_data <- replicate(
  10, better_simulate_si(
  params_inc, params_inf, params_iso, offset, max_shed, 4e4
)
  )

resized <- list()

for(i in 1:10){ # run multiple simulations to check variability
  
  unfiltered <- as.data.frame(sim_data[,i])
  
  filtered <- unfiltered[unfiltered$t_1 <= unfiltered$nu, ]
  
  # check simulated inf
  print(mean(unfiltered$t_1))
  print(sd(unfiltered$t_1))
  
  #print(mean(filtered[[i]]$nu))
  
  ## Make sure we have at least 500 rows
  idx <- sample(nrow(filtered), nsim_post_filter, replace = FALSE)
  resized[[i]] <- filtered[idx, ]
  
  # plot unfiltered, filtered, and filtered resampled to ensure 500 rows
  p[[i]] <- ggplot() +
    geom_density(
      aes(unfiltered$t_1), fill = "blue", col = NA, alpha = 0.2
    ) +
    geom_density(
      aes(resized[[i]]$t_1), fill = "red", col = NA, alpha = 0.2
    ) +
    geom_density(
      aes(filtered$t_1), fill = "red", col = NA, alpha = 0.2
    ) +
    geom_vline(xintercept = offset)+
    theme_minimal()
  
  max_si[i] <- ceiling(max(unfiltered$si))
}


for(i in 1:10){ # run multiple simulations to check variability
print(i)
y_vec <- seq(offset, max_si[i], by = 1)


###### grid likelihood
grid <- expand.grid(
  alpha1 = seq(1, 10, 0.5), beta1 = seq(1, 30, 0.5)
)

grid$normalised <- pmap_dbl(
  grid,
  function(alpha1, beta1) {
    out <- pmap_dbl(
      resized[[i]][, c("si", "nu")],
      function(si, nu) {
        full_model_lpdf(
          si, nu, max_shed, offset, 0, alpha1, beta1, alpha2, beta2,
          0.1, max_si[i], offset
          ) #-  normalising_constant(
          #y_vec, nu, max_shed, offset, 0, alpha1, beta1, alpha2, beta2,
          #0.1, max_si[i], offset
        #)
      }
    )
    sum(out)
  }
)

mle[[i]] <- grid[which(grid$normalised == max(grid$normalised)),]

mean[[i]] <- ((max_shed-offset)*(beta_shape1shape22muvar(mle[[i]]$alpha1, mle[[i]]$beta1)$mu))+offset

t_1_posterior[[i]] <- ((max_shed - offset) * rbeta(
  10000, shape1 = mle[[i]]$alpha1, shape2 = mle[[i]]$beta1))+offset

}

# plot the different posteriors onto one graph

t1_posterior_df <- data.frame(matrix(unlist(t_1_posterior), ncol=length(t_1_posterior), byrow=F))
t1_posterior_long <- reshape2::melt(t1_posterior_df)

# extract the different filtered t_1 data
filtered_t1 <- matrix(ncol = 10, nrow = length(resized[[1]]$t_1))

for(i in 1:10){
  
  filtered_t1[,i] <- resized[[i]]$t_1
}

# extract the different unfiltered t_1 data
unfiltered_t1 <- matrix(ncol = 10, nrow = length(sim_data[,1]$t_1))

for(i in 1:10){
  unfiltered_t1[,i] <- sim_data[,i]$t_1
}

filtered_t1_long <- reshape2::melt(as.data.frame(filtered_t1))
unfiltered_t1_long <- reshape2::melt(as.data.frame(unfiltered_t1))


q <- ggplot()+
  geom_density(data = unfiltered_t1_long,
    aes(value, group = variable), col = "blue", alpha = 0.2
  ) +
  geom_density(data = filtered_t1_long,
    aes(value, group = variable), colour = "red", alpha = 0.2
  ) +
  geom_density(data = t1_posterior_long,
    aes(value, group = variable), colour = "darkgreen"
    )+
  theme_minimal()
  


# print the means of each filtered dataset
mean_filtered <- vector(length = 10)
for(j in 1:10){
 mean_filtered[j] <- mean(resized[[j]]$t_1) 
}

# plot the different posteriors onto one graph

t1_posterior_df_n <- data.frame(matrix(unlist(t_1_posterior_norm), ncol=length(t_1_posterior_norm), byrow=F))
t1_posterior_long_n <- reshape2::melt(t1_posterior_df_n)

# extract the different filtered t_1 data
filtered_t1_n <- matrix(ncol = 10, nrow = length(resized[[1]]$t_1))

for(i in 1:10){
  
  filtered_t1_n[,i] <- resized[[i]]$t_1
}

# extract the different unfiltered t_1 data
unfiltered_t1_n <- matrix(ncol = 10, nrow = length(sim_data[,1]$t_1))

for(i in 1:10){
  unfiltered_t1_n[,i] <- sim_data[,i]$t_1
}

filtered_t1_long_n <- reshape2::melt(as.data.frame(filtered_t1_n))
unfiltered_t1_long_n <- reshape2::melt(as.data.frame(unfiltered_t1_n))


q_norm <- ggplot()+
  geom_density(data = unfiltered_t1_long_n,
               aes(value, group = variable), col = "blue", alpha = 0.2
  ) +
  geom_density(data = filtered_t1_long_n,
               aes(value, group = variable), colour = "red", alpha = 0.2
  ) +
  geom_density(data = t1_posterior_long_n,
               aes(value, group = variable), colour = "darkgreen"
  )+
  theme_minimal()





p3 <- ggplot(
  grid, aes(x = alpha1, beta1, fill = normalised)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd") +
  theme_minimal()



fits <- pmap(
  list(
    sim_data = with_recall_bias,
    param_inc = params_inc_all,
    param_offset = params_offset_all
  ),
  function(sim_data, param_inc, param_offset) {
    ## Choose a coarse grid here to make things faster
    si_vec <- seq(param_offset + 0.1 + 0.001, max_si, 1)
    fit <- stan(
      file = here::here("stan-models/full_model.stan"),
      data = list(
        N = length(sim_data$si),
        si = sim_data$si,
        nu = sim_data$nu,
        max_shed = max_shed,
        offset1 = param_offset,
        max_si = max(sim_data$si) + 0.001,
        min_si = param_offset, ## assuming the smallest incubation period is 0
        alpha2 = param_inc[["shape"]],
        beta2 = 1 / param_inc[["scale"]],
        width = width,
        M = length(si_vec),
        y_vec = si_vec,
        recall = 0
      ),
      chains = 2,
      iter = 3500,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
  }
)

## offset should be negative
posterior_inf_params <- function(fit, nsim, max_shed, offset) {
  params <- rstan::extract(fit)
  idx <- sample(length(params[[1]]), nsim, replace = TRUE)
  shape1 <- params[["alpha1"]][idx]
  shape2 <- params[["beta1"]][idx]
  out <- beta_shape1shape22muvar(shape1, shape2)
  out[["mu"]] <- max_shed * out[["mu"]] + offset
  out[["sigma2"]] <- max_shed * max_shed * out[["sigma2"]]
  out
}

out <- posterior_inf_params(fits[[1]], 5000, max_shed, -1)


params_inf_true <- params[[param_grid$params_inf[1]]]

p1 <- ggplot() +
  geom_density(aes(out[[1]]), fill = "red", col = NA, alpha = 0.3) +
  geom_vline(
    xintercept = params_inf_true$mean_inf
  ) +
  geom_vline(
    xintercept = c(
      mean(out[[1]]),
      mean(out[[1]]) - sd(out[[1]]),
      mean(out[[1]]) + sd(out[[1]])
    ), linetype = "dashed"
  ) +
  expand_limits(x = 0) +
  theme_minimal() +
  xlab("Mean infectious period")


p2 <- ggplot() +
  geom_density(aes(sqrt(out[[2]])), fill = "red", col = NA, alpha = 0.3) +
  geom_vline(
    xintercept = params_inf_true$sd_inf
  ) +
  geom_vline(
    xintercept = c(
      mean(sqrt(out[[2]])),
      mean(sqrt(out[[2]])) - sd(sqrt(out[[2]])),
      mean(sqrt(out[[2]])) + sd(sqrt(out[[2]]))
    ), linetype = "dashed"
  ) +
  expand_limits(x = 0) +
  theme_minimal() +
  xlab("SD infectious period")

p <- cowplot::plot_grid(p1, p2, ncol = 1)
