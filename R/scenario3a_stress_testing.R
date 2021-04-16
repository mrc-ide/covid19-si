prefix <- "3a_mix_stress_testing_sim_"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = "iso_par1",
  params_offset = c("offset1", "offset2"),
  params_pinvalid = c("pinvalid1", "pinvalid2"),
  stringsAsFactors = FALSE
)


index <- 1:nrow(param_grid)
param_grid <- param_grid[index, ]

params_inf_all <- pmap(
  list(
    params_inf = param_grid$params_inf,
    params_offset = param_grid$params_offset
  ),
  function(params_inf, params_offset) {
    out <- params_check[[params_inf]]
    offset <- params_check[[params_offset]]
    beta_muvar2shape1shape2(
      (out$mean_inf - offset) / (max_shed - offset),
      out$sd_inf^2 / (max_shed - offset)^2
    )
  }
)

params_inc_all <- map(
  param_grid$params_inc,
  function(params_inc) params_check[[params_inc]]
)

params_iso_all <- map(
  param_grid$params_iso,
  function(params_iso) params_check[[params_iso]]
)

params_offsets_all <- map(
  param_grid$params_offset,
  function(params_offset) params_check[[params_offset]]
)

params_pinv <- map(
  param_grid$params_pinv,
  function(params_pinv) params_check[[params_pinv]]
)

## unconditional
uncdtnl_data <- pmap(
  list(
    params_inf = params_inf_all,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offsets_all
  ),
  function(params_inf, params_inc, params_iso, params_offset) {
    ## ok to use nsim_post_filter here as we won't filter
    better_simulate_si(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      nsim_post_filter
    )
  }
)

## with -ve nu
simulated_data <- pmap(
  list(
   dat = uncdtnl_data,
   params_offset = params_offsets_all
  ),
  function(dat, params_offset) {
    toss <- runif(nrow(dat), 0, 1)
    idx <- toss < 0.02
    dat$nu[idx] <- runif(length(which(idx)), params_offset, 0)
    dat
  }
)

## conditional
## simulated_data <- map(
##   unconditional_data,
##   function(sim_data) {
##     sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]
##     ## Make sure we have at least 200 rows.
##     idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
##     sim_data[idx, ]
##   }
## )

invalid_si <- map(
  params_iso_all,
  function(params_iso) {
    invalid_si <- rbeta(
      nsim_post_filter, shape1 = alpha_invalid, shape2 = beta_invalid
    )

    invalid_iso <- rgamma(
      nsim_post_filter, shape = params_iso$shape, scale = params_iso$scale
    )
    data.frame(si = invalid_si, nu = invalid_iso)
  }
)


invalid_remapped <- map2(
  simulated_data,
  invalid_si,
  function(valid, invalid) {
    ## invalid SIs are draws from beta. Map them into
    ## min and max of valid SI
    ## Can make min_si here much smaller than offset
    ## to make it more like real data, and then the else statement
    ## in the model will take care of those very large negatives
    f <- map_into_interval(0, 1, min_invalid_si, max_invalid_si)
    invalid$si <- f(invalid$si)
    invalid
  }
)

mixed <- pmap(
  list(
    valid = simulated_data,
    invalid = invalid_remapped,
    params_pinv = param_grid$params_pinv
  ),
  function(valid, invalid, params_pinv) {
    pinvalid <- params_check[[params_pinv]]
    toss <- runif(nrow(valid))
    valid$type <- "valid"
    invalid$type <- "invalid"
    rbind(
      valid[toss > pinvalid , c("si", "nu", "type")],
      invalid[toss <= pinvalid ,c("si", "nu", "type")]
    )
  }
)

sampled <- pmap(
  list(sim_data = mixed,
       params_offset = params_offsets_all
    ),
  function(sim_data, params_offset) {
    ## Round for consistency with real data
    sim_data$si <- round(sim_data$si)
    sim_data$nu <- round(sim_data$nu)
    sim_data <- sim_data[sim_data$si > params_offset, ]
    sim_data <- sim_data[sim_data$nu > params_offset, ]
    idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
    sim_data[idx, ]
  }
)

outfiles <- glue::glue("data/{prefix}{seq_along(sampled)}_data.rds")
walk2(sampled, outfiles, function(x, y) saveRDS(x, y))

figs <- pmap(
  list(x = sampled, y = simulated_data, z = mixed, index = index),
  function(x, y, z, index){
    p <- ggplot() +
      geom_density(aes(x$si, fill = "red"), col = NA, alpha = 0.3) +
      geom_density(aes(y$si, fill = "blue"), col = NA, alpha = 0.3) +
      geom_density(aes(z$si, fill = "green"), col = NA, alpha = 0.3) +
      scale_fill_identity(
        breaks = c("red","blue", "green"),
        labels = c("Sampled", "Unconditional", "Mixed"),
        guide = "legend"
      ) +
      theme_minimal() +
      theme(legend.position = "top", legend.title = element_blank())
    ggsave(glue::glue("figures/{prefix}{index}_simulated.png"), p)
  }
)

iwalk(
  sampled, function(smpld, index) {
    p <- ggplot() +
      geom_histogram(aes(smpld$nu), binwidth = 1) +
      xlab("Delay from onset to isolation") +
      theme_minimal()
    cowplot::save_plot(glue::glue("figures/{prefix}{index}_nu.png"), p)
  }
)

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offsets_all,
    params_iso = params_iso_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, params_iso, sim_data, index) {
    si_vec <- seq(params_offset + 0.5, max_valid_si, 1)
    width <- 0.1
    sim_data <- arrange(sim_data, nu)
    fit_4a <- stan(
      file = here::here("stan-models/scenario3a_mixture_general.stan"),
      data = list(
        N = length(sim_data$si),
        si = sim_data$si,
        nu = sim_data$nu,
        max_shed = max_shed,
        offset1 = params_offset,
        alpha2 = params_inc[["shape"]],
        beta2 = 1 / params_inc[["scale"]],
        alpha_invalid = alpha_invalid,
        beta_invalid = beta_invalid,
        max_valid_si = max_valid_si,
        min_valid_si = params_offset,
        min_invalid_si = min_invalid_si,
        max_invalid_si = max_valid_si,
        width = width,
        M = length(si_vec),
        si_vec = si_vec,
        first_valid_nu = 1
      ),
      chains = 2, iter = 1000,
      seed = 42,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
    outfile <- glue::glue("stanfits/{prefix}_{index}.rds")
    saveRDS(fit_4a, outfile)
    fit_4a
  }
)
## index <- 1:nrow(param_grid)
## infiles <- glue::glue("stanfits/{prefix}_{index}.rds")
## fits <- map(infiles, readRDS)

process_fits <- pmap_dfr(
  list(
    fit = fits,
    offset = params_offsets_all
  ),
  function(fit, offset) {
    samples <- rstan::extract(fit)
    out <- mu_sd_posterior_distr(samples, max_shed, offset)
    pivot_wider(
      out, names_from = c("param", "var"), values_from = "val"
    )
  }, .id = "sim"
)

process_fits$true_mean <- map_dbl(
  params_check[param_grid$params_inf], function(x) x$mean_inf
)

process_fits$true_sd <- map_dbl(
  params_check[param_grid$params_inf], function(x) x$sd_inf
)

process_fits$incubation <- map_dbl(
  params_check[param_grid$params_inc],
  function(x) {
    epitrix::gamma_shapescale2mucv(x$shape, x$scale)[["mu"]]
  }
)

process_fits$isolation <- map_dbl(
  params_check[param_grid$params_iso],
  function(x) {
    epitrix::gamma_shapescale2mucv(x$shape, x$scale)[["mu"]]
  }
)

process_fits$pinvalid <- unlist(params_check[param_grid$params_pinv])
process_fits$offset <- unlist(params_check[param_grid$params_offset])

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

saveRDS(process_fits, glue::glue('stanfits/{prefix}processed_fits.rds'))

######### Posterior SI Distribution
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
    params_offset = params_offsets_all
  ),
  function(params_inf, params_inc, params_iso, params_offset) {
    sim_data <- better_simulate_si(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      nsim_pre_filter
    )
    sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]
    ##sim_data <- sim_data[abs(sim_data$si) > 0.1, ]
    ## Make sure we have at least 200 rows.
    idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
    sim_data[idx, ]
  }
)

mixed <- pmap(
  list(
    valid = posterior_si,
    invalid = invalid_remapped,
    pinvalid = process_fits$pinvalid_best
  ),
  function(valid, invalid, pinvalid) {
    toss <- runif(nrow(valid))
    valid$type <- "valid"
    invalid$type <- "invalid"
    rbind(
      valid[toss > pinvalid , c("si", "nu", "type")],
      invalid[toss <= pinvalid ,c("si", "nu", "type")]
    )
  }
)

sampled <- map(
  mixed, function(sim_data) {
    ## Round for consistency with real data
    sim_data$si <- round(sim_data$si)
    sim_data$nu <- round(sim_data$nu)
    sim_data <- sim_data[sim_data$si > 0, ]
    sim_data <- sim_data[sim_data$nu > 0, ]
    idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
    sim_data[idx, ]
  }
)

outfiles <- glue::glue("stanfits/{prefix}{seq_along(mixed)}_posterior.rds")
walk2(mixed, outfiles, saveRDS)
