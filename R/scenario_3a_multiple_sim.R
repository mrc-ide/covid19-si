## Set up params grid
prefix <- "3a_mix_sim_"
param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = c("iso_par1", "iso_par2"),
  params_offset = c("offset1", "offset2", "offset3"),
  params_pinvalid = c("pinvalid1", "pinvalid2", "pinvalid3"),
  stringsAsFactors = FALSE
)
param_grid <- param_grid[1, ]


params_inf_all <- pmap(
  list(
    params_inf = param_grid$params_inf,
    params_offset = param_grid$params_offset
  ),
  function(params_inf, params_offset) {
    out <- params[[params_inf]]
    ## The whole shifting in simulation will shift mu,
    ## so pass larger my to simulate function/
    offset <- params[[params_offset]]
    beta_muvar2shape1shape2(
    (out$mean_inf - offset) / (max_shed - offset),
    out$sd_inf^2 / (max_shed - offset)^2
    )
  }
)

params_inc_all <- map(
  param_grid$params_inc,
  function(params_inc) {
    out <- params[[params_inc]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_iso_all <- map(
  param_grid$params_iso,
  function(params_iso) {
    out <- params[[params_iso]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_offset_all <- map(
  param_grid$params_offset,
  function(params_offset) params[[params_offset]]
)

params_pinvalid_all <- map(
  param_grid$params_pinvalid,
  function(params_pinvalid) params[[params_pinvalid]]
)

simulated <- pmap(
  list(
    params_inc = params_inc_all,
    params_inf = params_inf_all,
    params_iso = params_iso_all,
    params_offset = params_offsets_all,
    params_pinvalid = params_pinvalid_all
  ),
  function(params_inc, params_inf, params_iso, params_offset,
           params_pinvalid) {
    min_si <- params_offset
    simulate_3a_mix(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_pre_filter, alpha_invalid, beta_invalid,
      min_si, max_si
    )
  }
)

sampled <- map(
  simulated, function(df) {
    ## might have NA from when pinvalid = 0
    out <- df[["simulated_si"]]
    out <- out[complete.cases(out), ]
    idx <- sample(nsim_pre_filter, nsim_post_filter, replace = TRUE)
    out[idx, ]
  }
)


fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offsets_all,
    sim_data = sampled,
    index = seq_along(params_inc_all)
  ),
  function(params_inc, params_offset, sim_data, index) {
    ## Rounding now to check things
    sim_data$si <- round(sim_data$si)
    fit_3a <- stan(
      file = here::here("stan-models/scenario3a_mixture_general.stan"),
      data = list(
        N = length(sim_data$si),
        si = sim_data$si,
        max_shed = max_shed,
        offset1 = params_offset,
        alpha2 = params_inc[["shape"]],
        beta2 = 1 / params_inc[["scale"]],
        alpha_invalid = alpha_invalid,
        beta_invalid = beta_invalid,
        max_si = max_si,
        min_si = params_offset + 0.001,
        width = width
      ),
      chains = 3,
      iter = 2000,
      verbose = TRUE
    )
    ## Save it now in case R crashes
    saveRDS(fit, glue::glue("stanfits/{prefix}{index}.rds"))
    fit_3a
  }
)



## fits <- map(infiles, readRDS)
posterior_mu_sd <- map(
  fits,
  function(fit) {
    if (nrow(as.data.frame(fit)) == 0) return(NULL)
    best_params <- map_estimates(fit)
    ## Same transformation on fitted params
    shifted_params <- beta_shape1shape22muvar(
      best_params[["alpha1"]], best_params[["beta1"]]
    )
    list(
      mu = (shifted_params$mu * (max_shed - offset)) + offset,
      sigma2 = (shifted_params$sigma2 * (max_shed - offset)^2)
    )
  }
)

posterior_sim <- pmap(
  list(
    fit = fits,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offsets_all
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

posterior_plots <- pmap(
  list(
    posterior_si = posterior_sim,
    sim_data = simulated_data
  ),
  function(posterior_si, sim_data) {
    if (is.null(posterior_si)) return(NULL)
    p <- ggplot() +
      geom_density(
        aes(sim_data$si, fill = "blue"), alpha = 0.3, col = NA
      ) +
      geom_density(
        aes(posterior_si$si, fill = "red"), alpha = 0.3, col = NA
      ) +
      scale_fill_identity(
        guide = "legend",
        labels = c("Simulated", "Posterior"),
        breaks = c("blue", "red")
      ) +
      xlab("Serial Interval") +
      ylab("Probability Density") +
      theme_minimal() +
      theme(legend.title = element_blank())
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
    sim_data = simulated_data
  ),
  function(posterior_si, sim_data) {
    if (is.null(posterior_si)) return(NULL)
    p <- ggplot() +
      geom_density(
        aes(sim_data$t_1, fill = "blue"), alpha = 0.3, col = NA
      ) +
      geom_density(
        aes(posterior_si$t_1, fill = "red"), alpha = 0.3, col = NA
      ) +
      scale_fill_identity(
        guide = "legend",
        labels = c("Simulated", "Posterior"),
        breaks = c("blue", "red")
      ) +
      xlab("Infectious Profile") +
      ylab("Probability Density") +
      theme_minimal() +
      theme(legend.title = element_blank())
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

ggsave("figures/scenario3a_si_multiple_sim.png", p)




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

ggsave("figures/scenario3a_t1_multiple_sim_params.png", p)

## alternative plot comparing the simulated and posterior infectious profile

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
  geom_point(
    aes(y = mu_simulated, x = sim), col = "red"
  ) +
  geom_linerange(
    aes(sim, ymin = ymin_simulated, ymax = ymax_simulated),
    size = 1, linetype = "dashed", col = "red",
    position = position_dodge(width = 0.4)
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
ggsave("figures/scenario3a_t1_multiple_sim_params.png", p)
