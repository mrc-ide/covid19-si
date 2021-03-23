prefix <- "4a_mix_with_stress_test_misspec_long_sim"

misspec_offset <- -7

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = "iso_par1",
  params_offset = c("offset1", "offset2"),
  params_pinvalid = c("pinvalid1", "pinvalid2"),
  params_beta = "beta1",
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
  param_grid$params_pinv, function(params_pinv) params_check[[params_pinv]]
)

params_beta <- map(
  param_grid$params_beta, function(params_beta) params_check[[params_beta]]
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
    better_simulate_si(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      nsim_pre_filter
    )
  }
)


## conditional
simulated_data <- map(
  uncdtnl_data,
  function(sim_data) {
    sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]
    ## Make sure we have at least 200 rows.
    idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
    sim_data[idx, ]
  }
)

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

## with -ve nu
mixed <- pmap(
  list(
   dat = mixed,
   params_offset = params_offsets_all
  ),
  function(dat, params_offset) {
    toss <- runif(nrow(dat), 0, 1)
    for(i in 1:(nrow(dat))){
      if(toss[i] < 0.02) {
        dat$nu[i] <- runif(1, params_offset, 0)
      }
    }
    dat[,]
  }
)

## Sample with recall
mixed <- pmap(
  list(
   dat = mixed,
   param_beta = params_beta
  ),
  function(dat, param_beta) {
    precall <- exp(-param_beta * abs(dat$si - dat$nu))
    idx <- sample(1:nrow(dat), nrow(dat), precall, replace = TRUE)
    dat[idx,]
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

outfiles <- glue::glue("data/{prefix}_{seq_along(mixed)}data.rds")
walk2(mixed, outfiles, function(x, y) saveRDS(x, y))

figs <- pmap(
  list(
    x = uncdtnl_data, y = simulated_data, z = sampled,
    index = index
  ),
  function(x, y, z, index){
    p <- ggplot() +
      geom_density(aes(x$si, fill = "red"), col = NA, alpha = 0.3) +
      geom_density(aes(y$si, fill = "blue"), col = NA, alpha = 0.3) +
      geom_density(aes(z$si, fill = "green"), col = NA, alpha = 0.3) +
      scale_fill_identity(
        breaks = c("red", "blue", "green"),
        labels = c("Conditional on nu", "Unconditional", "Mixed"),
        guide = "legend"
      ) +
      theme(legend.position = "top", legend.title = element_blank())
    ggsave(glue::glue("figures/{prefix}{index}_simulated.png"))

  }
)


fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = misspec_offset,
    params_iso = params_iso_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, params_iso, sim_data, index) {
    si_vec <- seq(params_offset + 0.5, max_valid_si, 1)
    width <- 0.1
    sim_data <- arrange(sim_data, nu)
    fit_4a <- stan(
      file = here::here("stan-models/scenario4a_mixture.stan"),
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
      chains = 2, iter = 2000,
      seed = 42,
      verbose = FALSE
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

## debugging
## checked that constant for full model with recall = 0 is the same as
## constant for s4.
## normalising_constant(y_vec, 5, 21, -2, 0, 9.52381, 15.47619, 9, 1 / 0.6666, 0.1, 21, -2)
## s4_normalising_constant(5, 21, -2, 9.52381, 15.47619, 9, 1 / 0.6666, 21, 0.1)

process_fits <- pmap_dfr(
  list(
    fit = fits,
    offset = misspec_offset),
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
    axis.title.y = element_blank()
  )

cowplot::save_plot(glue::glue("figures/{prefix}_inf_mu.png"), p)

psd <- ggplot(process_fits) +
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
    axis.title.y = element_blank()
  )

cowplot::save_plot(glue::glue("figures/{prefix}_inf_sd.png"), psd)

ppinv <- ggplot(process_fits) +
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
    axis.title.y = element_blank()
  ) +
  ylim(0, 0.5)

cowplot::save_plot(glue::glue("figures/{prefix}_inf_pinv.png"), ppinv)
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
    params_offset = misspec_offset
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

outfiles <- glue::glue("data/{prefix}_{seq_along(mixed)}data.rds")
training <- map(outfiles, readRDS)


compare_si <- map2_dfr(
  training, sampled,
  function(x, y) {
    x$si <- round(x$si)
    y$nu <- round(x$nu)
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

psi <- ggplot(process_fits) +
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
    axis.title.y = element_blank(), legend.title = element_blank()
  )
cowplot::save_plot(glue::glue("figures/{prefix}_inf_si.png"), psi)
