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
    params_offset = params_offset_all,
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
        min_si = params_offset - 0.001,
        width = width
      ),
      chains = 1,
      iter = 1000,
      verbose = TRUE
    )
    ## Save it now in case R crashes
    saveRDS(fit_3a, glue::glue("stanfits/{prefix}{index}.rds"))
    fit_3a
  }
)

process_fits <- pmap_dfr(
  list(fit = fits, offset = params_offset_all,
       true_vals = param_grid$params_inf,
       pinvalid = params_pinvalid_all),
  function(fit, offset, true_vals, pinvalid) {
    true_vals <- params[[true_vals]]
    true_val_df <- data.frame(
        true_val = c(true_vals$mean_inf, true_vals$sd_inf, pinvalid),
        param = c("mu", "sd", "pinvalid")
    )
    samples <- rstan::extract(fit)
    pinvalid_posterior <- quantile_as_df(samples[["pinvalid"]])
    pinvalid_posterior$param <- "pinvalid"
    out <- mu_sd_posterior_distr(samples, max_shed, offset)
    out <- rbind(out, pinvalid_posterior) %>%
      tidyr::spread(var, val)

    left_join(out, true_val_df)
  }, .id = "sim"
)

outfile <- glue::glue("data/{prefix}params_posterior_distr.rds")
saveRDS(process_fits, outfile)

process_fits$sim <- as.integer(process_fits$sim)

x <- split(process_fits, process_fits$param)

x[[3]] <- arrange(x[[3]], true_val)

ggplot(x[[3]]) +
  geom_linerange(aes(x = sim, ymin = `25%`, ymax = `75%`)) +
  geom_point(aes(sim, `50%`)) +
  geom_point(aes(sim, true_val), shape = 4) +
  facet_wrap(~true_val, scales = "free_y", ncol = 1) +
  expand_limits(y = 0) +
  theme_minimal()
