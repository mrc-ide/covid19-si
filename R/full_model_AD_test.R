max_shed <- 21
mean_inf <- 6
sd_inf <- 2
mean_inc <- 3
sd_inc <- 1
offset <- 0
iso <- c(rep(2, 100),9) # over represent 2s as they get filtered out harshly 

params_inf <- beta_muvar2shape1shape2(
  (mean_inf-offset)/(max_shed-offset), sd_inf^2 /(max_shed-offset)^2
)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = sd_inc/ mean_inc
)

nsim_post_filter <- 1000

sim_data <- better_simulate_si(
  params_inc, params_inf, params_iso, offset, max_shed, 1e4
)

iso_vec <- sample(iso, size = length(sim_data$nu), replace = TRUE)

sim_data$nu <- iso_vec
unfiltered <- sim_data
sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]

## Make sure we have at least 1000 rows.
idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
sim_data <- sim_data[idx, ]

ggplot() +
  geom_density(
    aes(unfiltered$t_1), fill = "blue", col = NA, alpha = 0.2
  ) +
  geom_density(
    aes(sim_data$t_1), fill = "red", col = NA, alpha = 0.2
  ) +
  geom_vline(xintercept = offset)+
  theme_minimal()

# check simulated inf
mean(unfiltered$t_1)
sd(unfiltered$t_1)

mean(sim_data$nu)



###### grid likelihood
grid <- expand.grid(
  alpha1 = seq(1, 10, 0.5), beta1 = seq(1, 20, 0.5)
)
## grid <- expand.grid(
##   alpha1 = params_inf$shape1, beta1 = seq(1, 30, 0.5)
## )

expose_stan_functions("stan-models/likelihoods.stan")

#max_si <- ceiling(max(sim_data$si))

max_si <- ceiling(max(unfiltered$si))

alpha2 <- params_inc$shape
beta2 <- 1 / params_inc$scale

y_vec <- seq(offset, max_si, by = 1)

grid$normalised <- pmap_dbl(
  grid,
  function(alpha1, beta1) {
    out <- pmap_dbl(
      sim_data[, c("si", "nu")],
      function(si, nu) {
        full_model_lpdf(
          si, nu, max_shed, offset, 0, alpha1, beta1, alpha2, beta2,
          0.1, max_si, offset
        )
        #-  normalising_constant(
        #y_vec, nu, max_shed, offset, 0, alpha1, beta1, alpha2, beta2,
        #0.1, max_si, offset
        #)
      }
    )
    sum(out)
  }
)

mle_6_2_0_norm <- grid[which(grid$normalised == max(grid$normalised)),]

mean_6_2_0_norm <- ((max_shed-offset)*(beta_shape1shape22muvar(mle_6_2_0_norm$alpha1, mle_6_2_0_norm$beta1)$mu))+offset

t_1_posterior_norm <- ((max_shed - offset) * rbeta(
  10000, shape1 = mle_6_2_0_norm$alpha1, shape2 = mle_6_2_0_norm$beta1))+offset

ggplot() +
  geom_density(
    aes(unfiltered$t_1), fill = "blue", col = NA, alpha = 0.2
  ) +
  geom_density(
    aes(sim_data$t_1), fill = "red", col = NA, alpha = 0.2
  ) +
  geom_density(
    aes(t_1_posterior_norm), fill = "green", col = NA, alpha = 0.2
  ) +
  geom_vline(xintercept = offset)+
  theme_minimal()


p2 <- ggplot(
  grid, aes(x = alpha1, beta1, fill = normalised)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd") +
  theme_minimal()
