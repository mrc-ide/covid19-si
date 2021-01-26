## multisimulation scenario 4 using full model: 5 simulated datasets 

#1 SIMULATE THE DATA

# true parameter values
max_shed <- 21
mean_inf <- 4
sd_inf <- 2
mean_inc <- 2
sd_inc <- 1
mean_iso <- 2
sd_iso <- 1
param_offset <- -1


params_inf <- beta_muvar2shape1shape2(
  (mean_inf-param_offset)/(max_shed-param_offset), sd_inf^2 /(max_shed-param_offset)^2
)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = sd_inc/ mean_inc
)

alpha2 <- params_inc$shape
beta2 <- 1 / params_inc$scale

params_iso <- epitrix::gamma_mucv2shapescale(
  mu = mean_iso, cv = sd_iso / mean_iso
)

nsim_post_filter <- 4000
n_datasets <- 3

sim_data <- replicate(
  n_datasets, better_simulate_si(
    params_inc, params_inf, params_iso, param_offset, max_shed, 2e5
  )
)

filtered_500 <- list()
p <- list()
max_si <- vector(length = n_datasets)

for(i in 1:n_datasets){ # run multiple simulations to check variability
  
  unfiltered <- as.data.frame(sim_data[,i])
  
  filtered <- unfiltered[unfiltered$t_1 <= unfiltered$nu, ]
  
  # check simulated inf
  print(mean(unfiltered$t_1))
  print(sd(unfiltered$t_1))
  
  #print(mean(filtered[[i]]$nu))
  
  ## Make sure we have at least 500 rows
  idx <- sample(nrow(filtered), nsim_post_filter, replace = FALSE)
  filtered_500[[i]] <- filtered[idx, ]
  
  max_si[i] <- ceiling(max(unfiltered$si))
}

#2 CHECK WHAT THE SIMULATED DATA LOOKS LIKE

# extract the different filtered t_1 data
filtered_t1 <- matrix(ncol = n_datasets, nrow = length(filtered_500[[1]]$t_1))

for(i in 1:n_datasets){
  
  filtered_t1[,i] <- filtered_500[[i]]$t_1
}

# extract the different unfiltered t_1 data
unfiltered_t1 <- matrix(ncol = n_datasets, nrow = length(sim_data[,1]$t_1))

for(i in 1:n_datasets){
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
  theme_minimal()

min_sim <- vector(length = n_datasets)
max_sim <- vector(length = n_datasets)
mean_sim <- vector(length = n_datasets)
median_sim <- vector(length = n_datasets)
sd_sim <- vector(length = n_datasets)

for(i in 1:n_datasets){
  min_sim[i] <- min(filtered_t1[,i])
  max_sim[i] <- max(filtered_t1[,i])
  mean_sim[i] <- mean(filtered_t1[,i])
  median_sim[i] <- median(filtered_t1[,i])
  sd_sim[i] <- sd(filtered_t1[,i])
}

#########################
## 3 FIT THE DATA SETS ##
#########################

fit <- list()


for(i in 1:n_datasets){
  
  si_vec <- seq(param_offset, max_si[i], by = 1)
  
  
  fit[[i]] <- stan(
  file = here::here("stan-models/full_model_no_norm.stan"),
  data = list(
    N = length(filtered_500[[i]]$si),
    si = filtered_500[[i]]$si,
    nu = filtered_500[[i]]$nu,
    max_shed = max_shed,
    offset1 = param_offset,
    max_si = max(sim_data[,i]$si) + 0.001,
    min_si = param_offset, ## assuming the smallest incubation period is 0
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    width = 0.1,
    #M = length(si_vec),
    #y_vec = si_vec,
    recall = 0
  ),
  chains = 2,
  iter = 4000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)
}

test_fit <- ggmcmc(ggs(fit[[3]]), here::here("figures/test3_4000.pdf"))

fitted_params <- list()
max_index <- list()
fitted_max <- list()
shape1_max <- vector(length = n_datasets)
shape2_max <- vector(length = n_datasets)
mean <- vector(length = n_datasets)
sd <- vector(length = n_datasets)
x <- matrix(nrow = 100000, ncol = n_datasets, byrow = FALSE)

for(i in 1:n_datasets){
 fitted_params[[i]] <- rstan::extract(fit[[i]]) 
 max_index[[i]] <- which(fitted_params[[i]]$lp__==max(fitted_params[[i]]$lp__))
 fitted_max[[i]] <- c(alpha1 = fitted_params[[i]]$alpha1[max_index[[i]]], 
                beta1 = fitted_params[[i]]$beta1[max_index[[i]]])
shape1_max[i] <- fitted_max[[i]]["alpha1"]
shape2_max[i] <- fitted_max[[i]]["beta1"]
mean[i] <- ((max_shed-offset)*beta_shape1shape22muvar(shape1_max[i], shape2_max[i])$mu)+offset
sd[i] <- (max_shed-offset)*sqrt(beta_shape1shape22muvar(shape1_max[i], shape2_max[i])$sigma2)
x[,i] <- (max_shed *
        rbeta(
          n = 100000,
          shape1 = shape1_max[i],
          shape2 = shape2_max[i]
        ))+offset
}

saveRDS(fitted_params, file = "fitted_params_check_4000.rds") 


# present the fits alongside the filtered and unfiltered data

# print the means of each filtered dataset
mean_filtered <- vector(length = n_datasets)
for(j in 1:n_datasets){
  mean_filtered[j] <- mean(filtered_500[[j]]$t_1) 
}

# plot the different posteriors onto one graph

t1_posterior_df <- data.frame(x)
t1_posterior_long <- reshape2::melt(t1_posterior_df)


q <- ggplot()+
  geom_density(data = unfiltered_t1_long,
               aes(value, group = variable), col = "blue", alpha = 0.2
  ) +
  geom_density(data = filtered_t1_long,
               aes(value, group = variable), colour = "red", alpha = 0.2
  ) +
  geom_density(data = t1_posterior_long,
               aes(value, group = variable), col = "forestgreen", alpha = 0.2)+
  theme_minimal()


ggsave("figures/recovered_inf_4000.png", q)

# plot the mean and sd
recov <- data.frame(mean, sd, tag = "recov")
true <- data.frame(mean = mean_inf, sd = sd_inf, tag = "true")
df <- rbind(recov, true)
df <- reshape2::melt(df, id.vars = "tag")

ggplot()+
  geom_point(data = df, aes(x = variable, y = value, col = tag, pch = tag), size = 3)+
  theme_minimal()

# what about the filtered data is leading to the variance?

tabl <- data.frame(mean, sd, mean_sim, median_sim, min_sim, max_sim, sd_sim)
