## Set up params grid
prefix <- "3a_mix_with_normalisation_sim_"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = "iso_par1", ## S3 doesn't care about Isolation
  params_offset = c("offset1", "offset2", "offset3"),
  params_pinvalid = c("pinvalid1", "pinvalid2", "pinvalid3"),
  stringsAsFactors = FALSE
)

## Process stuck after running 44th row forever.
## Restarting
##index <- 45:nrow(param_grid)
index <- 1:nrow(param_grid)
param_grid <- param_grid[index, ]
## param_grid <- tail(param_grid, 1)

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
    params_offset = params_offset_all,
    params_pinvalid = params_pinvalid_all
  ),
  function(params_inc, params_inf, params_iso, params_offset,
           params_pinvalid) {
    min_si <- params_offset
    simulate_3a_mix(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_pre_filter, alpha_invalid, beta_invalid,
      min_invalid_si, max_invalid_si
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

outfiles <- glue::glue("data/{prefix}_{seq_along(sampled)}data.rds")
walk2(sampled, outfiles, function(x, y) saveRDS(x, y))



fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offset_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, sim_data, index) {
    ## Rounding now to check things
    si_vec <- seq(params_offset + 0.5, max_valid_si, 1)
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
        max_valid_si = max_valid_si,
        min_valid_si = params_offset,
        max_invalid_si = max_invalid_si,
        min_invalid_si = min_invalid_si,
        width = width,
        M = length(si_vec),
        si_vec = si_vec
      ),
      seed = 42,
      verbose = TRUE
    )
    ## Save it now in case R crashes
    saveRDS(fit_3a, glue::glue("stanfits/{prefix}{index}.rds"))
    fit_3a
  }
)


## index <- 1:nrow(param_grid)
## infiles <- glue::glue("stanfits/{prefix}{index}.rds")
## fits <- map(infiles, readRDS)


process_fits <- pmap_dfr(
  list(
    fit = fits,
    offset = params_offset_all),
  function(fit, offset) {
    samples <- rstan::extract(fit)
    out <- mu_sd_posterior_distr(samples, max_shed, offset)
    pivot_wider(
      out, names_from = c("param", "var"), values_from = "val"
    )
  }, .id = "sim"
)

process_fits$true_mean <- map_dbl(
  params[param_grid$params_inf], function(x) x$mean_inf
)

process_fits$true_sd <- map_dbl(
  params[param_grid$params_inf], function(x) x$sd_inf
)

process_fits$incubation <- map_dbl(
  params[param_grid$params_inc], function(x) x$mean_inc
)

process_fits$isolation <- map_dbl(
  params[param_grid$params_iso], function(x) x$mean_iso
)

process_fits$pinvalid <- unlist(params[param_grid$params_pinv])
process_fits$offset <- unlist(params[param_grid$params_offset])

est_pinvalid <- map_dfr(
  fits,
  function(fit) {
    samples <- rstan::extract(fit)
    out <- quantile_as_df(samples[["pinvalid"]])
    idx <- which.max(samples[["lp__"]])
    best <- samples[["pinvalid"]][idx]
    out <- rbind(out, data.frame(var = "best", val = best))
    out$param <- "pinvalid"
    pivot_wider(
      out, names_from = c("param", "var"), values_from = "val"
    )
  }, .id = "sim"
)

process_fits <- left_join(process_fits, est_pinvalid, by = "sim")


process_fits <- mutate_if(process_fits, is.numeric, ~ round(., 2))

p <- ggplot(process_fits) +
  geom_point(aes(sim, `mu_50%`)) +
  geom_linerange(aes(x = sim, ymin = `mu_25%`, ymax = `mu_75%`)) +
  geom_point(aes(sim, true_mean), shape = 4) +
  facet_grid(
    pinvalid ~ offset, scales = "free", labeller = label_both
  ) +
  ylab("Mean infectious period") +
  theme_minimal() +
  theme(
    legend.position = "top", axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
ggsave(glue::glue("figures/{prefix}mean_inf.png"), p)

p <- ggplot(process_fits) +
  geom_point(aes(sim, `sd_50%`)) +
  geom_linerange(aes(x = sim, ymin = `sd_25%`, ymax = `sd_75%`)) +
  geom_point(aes(sim, true_sd), shape = 4) +
  facet_grid(
    pinvalid ~ offset, scales = "free", labeller = label_both
  ) +
  ylab("SD infectious period") +
  theme_minimal() +
  theme(
    legend.position = "top", axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
ggsave(glue::glue("figures/{prefix}sd_inf.png"), p)

p <- ggplot(process_fits) +
  geom_point(aes(sim, `pinvalid_50%`)) +
  geom_linerange(aes(x = sim, ymin = `pinvalid_25%`, ymax = `pinvalid_75%`)) +
  geom_point(aes(sim, pinvalid), shape = 4) +
  facet_grid(
    pinvalid ~ offset, scales = "free", labeller = label_both
  ) +
  ylab("pinvalid") +
  theme_minimal() +
  theme(
    legend.position = "top", axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylim(0, 1)

ggsave(glue::glue("figures/{prefix}pinvalid.png"), p)


######### Posterior SI Distribution
training <- sampled
params_inf_post <- map(
  fits, function(fit) {
    samples <- rstan::extract(fit)
    idx <- which.max(samples[["lp__"]])
    list(
      shape1 = samples[["alpha1"]][idx],
      shape2 = samples[["beta1"]][idx]
    )
  }
)

posterior_si <- pmap(
  list(
    params_inf = params_inf_post,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offset_all,
    params_pinvalid = params_pinvalid_all
  ),
  function(params_inf, params_inc, params_iso, params_offset,
           params_pinvalid) {
    simulate_3a_mix(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_pre_filter, alpha_invalid, beta_invalid,
      min_invalid_si, max_si
    )
  }
)


sampled <- map(
  posterior_si, function(df) {
    out <- df[["simulated_si"]]
    out <- out[complete.cases(out), ]
    idx <- sample(nsim_pre_filter, nsim_post_filter, replace = TRUE)
    out[idx, ]
  }
)





compare_si <- map2_dfr(
  training, sampled,
  function(x, y) {
    x$si <- round(x$si)
    trng_si <- quantile_as_df(x$si)
    trng_si$param <- "training"
    post_si <- quantile_as_df(y$si)
    post_si$param <- "posterior"
    rbind(trng_si, post_si) %>%
      pivot_wider(
      names_from = c("param", "var"), values_from = "val"
    )
  }, .id = "sim"
)

process_fits <- left_join(process_fits, compare_si, by = "sim")
process_fits <- arrange(process_fits, `training_50%`)
process_fits$sim <- factor(process_fits$sim, process_fits$sim)

p <- ggplot(process_fits) +
  geom_point(aes(sim, `training_50%`, col = "red")) +
  geom_linerange(
    aes(
      x = sim, ymin = `training_25%`, ymax = `training_75%`,
      col = "red")
  ) +
  geom_point(
    aes(sim, `posterior_50%`, col = "blue"),
    position = position_nudge(x = 0.5)
  ) +
  geom_linerange(
    aes(
      x = sim, ymin = `posterior_25%`, ymax = `posterior_75%`,
      col = "blue"), position = position_nudge(x = 0.5)
  ) +
  scale_color_identity(
    breaks = c("red", "blue"), labels = c("Training", "Posterior"),
    guide = "legend"
  ) +
  facet_grid(
    pinvalid ~ offset, scales = "free", labeller = label_both
  ) +
  ylab("Serial Interval") +
  theme_minimal() +
  theme(
    legend.position = "top", axis.text.x = element_blank(),
    axis.title.x = element_blank(), legend.title = element_blank()
  )
ggsave(glue::glue("figures/{prefix}si.png"), p)
