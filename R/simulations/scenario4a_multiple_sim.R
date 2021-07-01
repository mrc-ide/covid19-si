max_shed <- 21
mean_inf <- 4
sd_inf <- 1
mean_inc <- 2
sd_inc <- 1
## very short isolation
mean_iso <- 5
sd_iso <- 2
offset <- -3

params_inf <- beta_muvar2shape1shape2(
  (mean_inf - offset) / (max_shed - offset),
  sd_inf^2 /(max_shed - offset)^2
)

params_inc <- epitrix::gamma_mucv2shapescale(
  mu = mean_inc, cv = sd_inc/ mean_inc
)

alpha2 <- params_inc$shape
beta2 <- 1 / params_inc$scale

params_offset_all <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, params_offset) {
    params[[params_offset]]
  }
)

nsim_post_filter <- 500
simulated_data <- pmap(
  list(
    params_inf = params_inf_all,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offset_all
  ),
  function(params_inf, params_inc, params_iso, params_offset) {
    sim_data <- better_simulate_si(
      params_inc, params_inf, params_iso, params_offset, max_shed, nsim
    )
    sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]
    sim_data <- sim_data[abs(sim_data$si) > 0.1, ]
    ## Make sure we have at least 200 rows.
    idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
    sim_data[idx, ]
  }
)

filtered_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]

idx <- sample(nrow(filtered_data), nsim_post_filter)
filtered_data <- filtered_data[idx, ]

prefix <- "4a_sim_"

p <- ggplot() +
  geom_density(
    aes(filtered_data$t_1, fill = "red"), alpha = 0.3, col = NA
  ) +
  geom_density(
    aes(sim_data$t_1, fill = "blue"), alpha = 0.3, col = NA
  ) +
  scale_fill_identity(
    breaks = c("red", "blue"),
    labels = c("secondary infection before isolation", "All"),
    guide = "legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank()) +
  xlab("Infectious Period")


## Without normalisation
width <- 0.1
max_si <- ceiling(max(sim_data$si)) + 1
si_vec <- seq(offset + 0.1 + 0.001, max_si, 1)

fit <- stan(
  ##file = here::here("stan-models/scenario4a.stan"),
  file = here::here("stan-models/scenario4a_with_normalisation.stan"),
  data = list(
    N = length(filtered_data$si),
    si = filtered_data$si,
    nu = filtered_data$nu,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    offset1 = offset,
    width = width,
    M = length(si_vec),
    y_vec = si_vec
  ),
  chains = 3,
  iter = 2000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)

best_params <- map_estimates(fit)
shifted_params <- beta_shape1shape22muvar(
  best_params[["alpha1"]], best_params[["beta1"]]
)

posterior_mu_sd <- list(
  mu = (shifted_params$mu * (max_shed - offset)) + offset,
  sigma2 = (shifted_params$sigma2 * (max_shed - offset)^2)
)


posterior_si <- better_simulate_si(
  params_inc,
  list(shape1 = best_params[["alpha1"]],
       shape2 = best_params[["beta1"]]),
  params_iso, offset, max_shed, nsim = 10000
)

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offset_all,
    sim_data = simulated_data
  ),
  function(params_inc, params_offset, sim_data) {

p <- ggplot() +
  geom_density(
    aes(sim_data$t_1, fill = "blue"), alpha = 0.3, col = NA
  ) +
  geom_density(
    aes(filtered_data$t_1, fill = "green"), alpha = 0.3, col = NA
  ) +
  geom_density(
     aes(posterior_si$t_1, fill = "red"), alpha = 0.3, col = NA
   ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Simulated", "Secondary infection before isolation",
               "Posterior"),
    breaks = c("blue", "green", "red")
  ) +
  xlab("Infectious Profile") +
  ylab("Probability Density") +
  theme_minimal() +
  theme(
    legend.title = element_blank(), legend.position = "top"
  )

ggsave("figures/posterior_t1_4a_sim_1.png", p)

posterior_samples <- rstan::extract(fit)
posterior_samples <- map(
  posterior_samples, function(x) x[sample(length(x), 1000)]
)

posterior_mu_distr <- beta_shape1shape22muvar(
  posterior_samples[["alpha1"]], posterior_samples[["beta1"]]
)

posterior_mu_distr <- list(
  mu = (posterior_mu_distr$mu * (max_shed - offset)) + offset,
  sigma2 = (posterior_mu_distr$sigma2 * (max_shed - offset)^2)
)


p_mu <- ggplot() +
  geom_density(
    aes(posterior_mu_distr$mu), fill = "red", col = NA, alpha = 0.3 ) +
  geom_vline(xintercept = mean_inf) +
  theme_minimal() +
  xlab("Mean Infectious Profile")

posterior_sim <- pmap(
  list(
    fit = fits,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offset_all
  ),
  function(fit, params_inc, params_iso, params_offset) {
    if (nrow(as.data.frame(fit)) == 0) return(NULL)
    best_params <- map_estimates(fit)
    better_simulate_si(
      params_inc,
      list(shape1 = best_params[["alpha1"]],
           shape2 = best_params[["beta1"]]),
      params_iso, params_offset, max_shed, nsim = 10000
    )
  }
)

p <- cowplot::plot_grid(p_mu, p_sd, ncol = 1)

ggsave("figures/posterior_params_4a_with_normalisation.png", p)



posterior_plots <- pmap(
  list(
    posterior_si = posterior_sim,
    sim_data = simulated_data,
    filtered = filtered_data
  ),
  function(posterior_si, sim_data, filtered) {
    if (is.null(posterior_si)) return(NULL)
    p <- ggplot() +
      geom_density(
        aes(sim_data$si, fill = "blue"), alpha = 0.3, col = NA
      ) +
      geom_density(
        aes(filtered$si, fill = "green"), alpha = 0.3, col = NA
      ) +
      geom_density(
        aes(posterior_si$si, fill = "red"), alpha = 0.3, col = NA
      ) +
      scale_fill_identity(
        guide = "legend",
        labels = c("Simulated All", "Simulated Filtered","Posterior"),
        breaks = c("blue", "green",  "red")
      ) +
      xlab("Serial Interval") +
      ylab("Probability Density") +
      theme_minimal() +
      theme(legend.title = element_blank(), legend.position = "top")
    p

  }
)

iwalk(
  posterior_plots,
  function(p, index) {
    ggsave(
    glue::glue("figures/posterior_si_{prefix}{index}.png"), p)
  }
)

posterior_t1_plots <- pmap(
  list(
    posterior_si = posterior_sim,
    sim_data = simulated_data,
    filtered = filtered_data
  ),
  function(posterior_si, sim_data, filtered) {
    if (is.null(posterior_si)) return(NULL)

    p

  }
)

iwalk(
  posterior_t1_plots,
  function(p, index) {
    ggsave(
    glue::glue("figures/posterior_t1_{prefix}{index}.png"), p)
  }
)


params_compare <- pmap_dfr(
  list(
    fit = fits,
    params_inf = param_grid$params_inf,
    params_offset = params_offset_all
  ),
  function(fit, params_inf, params_offset) {
    out <- params[[params_inf]]
    true_values <- data.frame(
      param = c("mu", "sd", "offset"),
      var = "simulated",
      val = c(out$mean_inf, out$sd_inf, params_offset)
    )
    if (nrow(as.data.frame(fit)) == 0) return(NULL)
    best_params <- map_estimates(fit)

    out <- rbeta(10000, best_params$alpha1, best_params$beta1)
    f <- map_into_interval(0, 1, params_offset, max_shed)
    out <- f(out)
    out <- data.frame(
      param = c("mu", "sd", "offset"),
      var = "posterior",
      val = c(mean(out), sd(out), params_offset)
    )
    out <- rbind(out, true_values)
    tidyr::pivot_wider(out, names_from = "param", values_from = "val")
  },
  .id = "sim"
)

si_compare <- pmap_dfr(
  list(
    posterior_si = posterior_sim,
    sim_data = simulated_data
  ),
  function(posterior_si, sim_data) {

    if (is.null(posterior_si)) return(NULL)
    simulated_summary <- quantile_as_df(sim_data$si)
    simulated_summary$category <- "Simulated"
    posterior_summary <- quantile_as_df(posterior_si$si)
    posterior_summary$category <- "Posterior"
    rbind(simulated_summary, posterior_summary)
  },
  .id = "sim"
)



si_compare <- tidyr::spread(si_compare, var, val)
si_compare$sim <- factor(
  si_compare$sim, levels = as.character(1:nrow(param_grid)), ordered = TRUE
)

p <- ggplot(si_compare) +
  geom_point(
    aes(sim, `50%`, col = category),
    size = 2,
    position = position_dodge(width = 0.3)
  ) +
  geom_linerange(
    aes(sim, ymin = `2.5%`, ymax = `97.5%`, col = category),
    size = 1.1,
    position = position_dodge(width = 0.3)
  ) +
  scale_x_discrete(
    breaks = as.character(1:nrow(param_grid))
  ) +
  scale_color_manual(
    values = c(Simulated = "blue", Posterior = "red"),
    labels = c("Simulated value", "Posterior Estimate")
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank()) +
  xlab("Simulation") +
  ylab("Serial Interval (Median and 95% CrI)")

ggsave("figures/scenario4a_si_multiple_sim.png", p)




params_compare$sim <- factor(
  params_compare$sim, levels = as.character(1:nrow(param_grid)), ordered = TRUE
)

params_compare$ymin <- params_compare$mu - params_compare$sd
params_compare$ymax <- params_compare$mu + params_compare$sd

x <- tidyr::pivot_wider(params_compare, names_from = "var",
                        values_from = c("mu", "sd", "ymin", "ymax"))

p <- ggplot(x) +
  geom_point(
    aes(sim, mu_posterior-offset),
    size = 2,
    position = position_dodge(width = 0.4)
  ) +
  geom_linerange(
    aes(sim, ymin = ymin_posterior-offset, ymax = ymax_posterior-offset),
    size = 1.1,
    position = position_dodge(width = 0.4)
  ) +
  geom_hline(
    aes(yintercept = mu_simulated), linetype = "dashed"
  ) +
  geom_hline(
    aes(yintercept = ymin_simulated), linetype = "dashed"
  ) +
  geom_hline(
    aes(yintercept = ymax_simulated), linetype = "dashed"
  ) +
  scale_x_discrete(
    breaks = as.character(1:nrow(param_grid))
  ) +
  ## scale_color_manual(
  ##   values = c(true = "blue", posterior = "red"),
  ##   labels = c("True value", "Posterior Estimate")
  ## ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank()) +
  xlab("Simulation") +
  ylab("Infectious profile")
  ##facet_wrap(~group, ncol = 2, scales = "free_x")

ggsave("figures/scenario4a_t1_multiple_sim_params.png", p)
