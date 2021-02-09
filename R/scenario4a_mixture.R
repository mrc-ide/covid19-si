prefix <- "4a_mix_with_normalisation_sim"

param_grid <- expand.grid(
  params_inf = c("inf_par2"),
  params_inc = c("inc_par2"),
  params_iso = c("iso_par1"),
  params_pinv = c("pinvalid1"),
  params_offset = c("offset3"),
  stringsAsFactors = FALSE
)


index <- 1:nrow(param_grid)
param_grid <- param_grid[index, ]

params_inf_all <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, params_pinv, params_offset) {
    out <- params[[params_inf]]
    beta_muvar2shape1shape2(
      out$mean_inf/max_shed, out$sd_inf^2 / max_shed^2
    )
  }
)

params_inc_all <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, params_pinv, params_offset) {
    out <- params[[params_inc]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_iso_all <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, params_pinv, params_offset) {
    out <- params[[params_iso]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_offsets_all <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, params_pinv, params_offset) {
    params[[params_offset]]
  }
)

params_pinv <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, params_pinv, params_offset) {
    params[[params_pinv]]
  }
)

simulated_data <- pmap(
  list(
    params_inf = params_inf_all,
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
    max_si <- max(valid$si)
    ## invalid SIs are draws from beta. Map them into
    ## min and max of valid SI
    ## Can make min_si here much smaller than offset
    ## to make it more like real data, and then the else statement
    ## in the model will take care of those very large negatives
    f <- map_into_interval(0, 1, min_invalid_si, max_si)
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
    pinvalid <- params[[params_pinv]]
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
walk2(mixed, outfiles, function(x, y) saveRDS(x, y))

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offsets_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, sim_data, index) {
    max_si <- max(sim_data$si) + 1
    width <- 0.1
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
        max_valid_si = max_si,
        min_valid_si = params_offset,
        min_invalid_si = min_invalid_si,
        max_invalid_si = max_si,
        width = width
      ),
      chains = 1, iter = 800,
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

## debugging
sim_data <- sampled[[1]]
si_vec <- seq(-2, 21, 1) + 1
nus <- sort(unique(sim_data$nu))
n_constant <- matrix(
  NA, nrow = length(si_vec), ncol = length(unique(sim_data$nu))
)

for (row in seq_along(si_vec)) {
  for (col in seq_along(nus)) {
    n_constant[row, col] <- exp(
      scenario4a_lpdf(
        si_vec[row], nus[col], max_shed = max_shed, offset1 = -2,
        alpha1 = 9.52381, beta1 = 15.47619, alpha2 = 9, beta2 = 1 / 0.67,
        width = 0.1
      )
    )
  }
}

normalise_valid <- colSums(n_constant)
normalise_total <- normalise_valid
##normalise_total <- log(normalise_valid)
names(normalise_total) <- nus
grid <- expand.grid(shape1 = seq(1, 15), shape2 = seq(1, 20))
## This works, ll max at shape1 = 9, shape2 = 15
## True values 9.5, 15.5
ll <- pmap_dbl(
  grid,
  function(shape1, shape2) {
    out <- pmap_dbl(
      sim_data[, c("si", "nu")], function(si, nu) {
        scenario4a_lpdf(
          si, nu, max_shed, -2, shape1, shape2, alpha2 = 9,
          beta2 = 1 / 0.67, width = 0.1
        ) - normalise_total[[as.character(nu)]]
      }
    )
    sum(out)
  }
)

grid$ll <- ll

process_fits <- pmap_dfr(
  list(
    fit = fits,
    offset = params_offsets_all),
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

ggplot(process_fits) +
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


ggplot(process_fits) +
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

ggplot(process_fits) +
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
  ylim(0, 1)


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

ggplot(process_fits) +
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
