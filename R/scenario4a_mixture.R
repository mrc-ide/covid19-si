param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2", "inf_par3"),
  params_inc = c("inc_par1", "inc_par2", "inc_par2"),
  params_iso = c("iso_par1", "iso_par2", "iso_par3"),
  params_pinv = c("pinvalid1", "pinvalid2", "pinvalid3"),
  params_offset = "offset1",
  stringsAsFactors = FALSE
)
prefix <- "4a_mix_"
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

outfiles <- glue::glue("data/{prefix}_{seq_along(mixed)}data.rds")
walk2(mixed, outfiles, function(x, y) saveRDS(x, y))

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offsets_all,
    sim_data = mixed,
    index = seq_along(mixed)
  ),
  function(params_inc, params_offset, sim_data, index) {

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
      iter = 3000,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
    outfile <- glue::glue("stanfits/{prefix}_{index}.rds")
    saveRDS(fit_4a, outfile)
    fit_4a
  }
)

params_posterior <- map(
  fits, function(fit) {
    if (nrow(as.data.frame(fit)) == 0) return(NULL)
    else map_estimates(fit)
  }
)

simulated_posterior <- pmap(
  list(
    params_inf = params_posterior,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offsets_all
  ),
  function(params_inf, params_inc, params_iso, params_offset) {
    if (is.null(params_inf)) return(NULL)

    params_inf$shape1 <- params_inf$alpha1
    params_inf$shape2 <- params_inf$beta1
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

invalid_posterior <- map2(
  simulated_posterior,
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

mixed_posterior <- pmap(
  list(
    valid = simulated_posterior,
    invalid = invalid_posterior,
    params_pinv = params_posterior
  ),
  function(valid, invalid, params_pinv) {
    if (is.null(params_pinv)) return(NULL)

    pinvalid <- params_pinv$pinvalid
    toss <- runif(nrow(valid))
    valid$type <- "valid"
    invalid$type <- "invalid"
    rbind(
      valid[toss > pinvalid , c("si", "nu", "type")],
      invalid[toss <= pinvalid ,c("si", "nu", "type")]
    )
  }
)

si_compare <- map2_dfr(
  mixed,
  mixed_posterior,
  function(simulated, posterior) {
    x <- quantile_as_df(simulated$si)
    x$type <- "Simulated"
    y <- quantile_as_df(posterior$si)
    y$type <- "Posterior"
    rbind(x, y)
  }, .id = "sim"
)

x <- na.omit(si_compare)
x <- tidyr::spread(x, var, val)

x$sim <- factor(
  x$sim, levels = as.character(1:nrow(param_grid)),
  ordered = TRUE
)

ggplot(x) +
  geom_point(aes(sim, `50%`, col = type)) +
  geom_linerange(
    aes(sim, ymin = `2.5%`, ymax = `97.5%`, col = type)
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90)
  ) +
  ylab("SI Median and 95% CrI") +
  xlab("Simulation")

params_compare <- pmap_dfr(
  list(
    simulated = param_grid$params_inf,
    posterior = params_posterior,
    offset = param_grid$params_offset
  ),
  function(simulated, posterior, offset) {
    if (is.null(params_posterior)) return(NULL)
    simulated <- params[[simulated]]
    offset <- params[[offset]]
    shape1 <- posterior$alpha1
    shape2 <- posterior$beta1
    out <- beta_shape1shape22muvar(shape1, shape2)
    out$mu <- (max_shed * out$mu) - offset
    df <- data.frame(
      type = c("simulated", "simulated",
               "posterior", "posterior"),
      var = c("mu", "sd", "mu", "sd"),
      val = c(simulated[[1]], simulated[[2]],
              out$mu, max_shed * sqrt(out$sigma2))
    )
    df$offset <- offset
    tidyr::pivot_wider(df, names_from = c("type", "var"),
                       values_from = "val")
  }, .id = "sim"
)

params_compare$sim <- factor(
  params_compare$sim, levels = as.character(1:nrow(param_grid)),
  ordered = TRUE
)

ggplot(params_compare) +
  geom_point(aes(sim, posterior_mu)) +
geom_linerange(
  aes(sim, ymin = posterior_mu - posterior_sd,
      ymax = posterior_mu + posterior_sd)
) +
  geom_point(aes(sim, simulated_mu), col = "red") +
  facet_wrap(~offset) +
  theme_minimal() +
    theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90)
  )
