max_shed <- 21
mean_inf <- 6
sd_inf <- 2
mean_inc <- 3
sd_inc <- 1
## very short isolation
mean_iso <- 2
sd_iso <- 0.1
offset <- 0

params_inf <- beta_muvar2shape1shape2(
  (mean_inf-offset)/(max_shed-offset), sd_inf^2 /(max_shed-offset)^2
)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = sd_inc/ mean_inc
)

params_iso <- epitrix::gamma_mucv2shapescale(
  mu = mean_iso, cv = sd_iso / mean_iso
)

nsim_post_filter <- 1000

sim_data <- better_simulate_si(
  params_inc, params_inf, params_iso, offset, max_shed, 1e4
)

unfiltered <- sim_data
sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]


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

## Make sure we have at least 200 rows.
idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
sim_data <- sim_data[idx, ]

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
