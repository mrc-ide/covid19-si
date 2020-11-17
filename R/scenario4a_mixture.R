param_grid <- expand.grid(
  params_inf = "inf_par1",
  params_inc = "inc_par3",
  params_iso = "iso_par1",
  params_pinv = "pinvalid3",
  params_offset = "offset1",
  stringsAsFactors = FALSE
)

nsim <- 20000

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

nsim_post_filter <- 500
simulated_data <- pmap(
  list(
    params_inf = params_inf_all,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offsets_all
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
    ##
    f <- map_into_interval(0, 1, 0.5 * min(valid$si), 2 * max(valid$si))
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

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offsets_all,
    sim_data = mixed
  ),
  function(params_inc, params_offset, sim_data) {

    width <- 0.05
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
        max_si = max(sim_data$si) + 0.001,
        min_si = min(sim_data$si) - 0.001,
        width = width
      ),
      chains = 3,
      iter = 2000,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
    fit_4a
  }
)

