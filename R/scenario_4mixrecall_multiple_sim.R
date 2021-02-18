## Set up params grid
prefix <- "4a_recall_sim_no_norm"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = "inc_par2",
  params_iso = c("iso_par1", "iso_par2"),
  params_offset = "offset3", # just try with offset = -3 initially
  params_pinvalid = "pinvalid1", #just try with pinvalid = 0 initially
  params_recall = c("recall1", "recall2"),
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

params_recall_all <- map(
  param_grid$params_recall,
  function(params_recall) params[[params_recall]]
)

simulated <- pmap(
  list(
    params_inc = params_inc_all,
    params_inf = params_inf_all,
    params_iso = params_iso_all,
    params_offset = params_offset_all,
    params_pinvalid = params_pinvalid_all,
    params_recall = params_recall_all
  ),
  function(params_inc, params_inf, params_iso, params_offset,
           params_pinvalid, params_recall) {
    min_si <- params_offset
    simulate_4a_mix(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_post_filter, alpha_invalid, beta_invalid,
      min_si, max_si, params_recall
    )
  }
)

sampled <- map(
  simulated, function(df) {
    ## might have NA from when pinvalid = 0
    out <- df[["simulated_si_recalled"]]
    out <- out[complete.cases(out), ]
    #idx <- sample(nsim_pre_filter, nsim_post_filter, replace = TRUE) - I sampled in my sim function
    out
  }
)

#si_vec <- seq(-3, 40, by = 1)

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offset_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, sim_data, index) {
    ## Rounding now to check things
    sim_data$si <- round(sim_data$si)
    fit_4a_recall <- stan(
      file = here::here("stan-models/full_model_no_norm.stan"),
      data = list(
        N = length(sim_data$si),
        si = sim_data$si,
        nu = sim_data$nu,
        max_shed = max_shed,
        offset1 = params_offset,
        max_si = max_si,
        min_si = params_offset - 0.001,  
        alpha2 = params_inc[["shape"]],
        beta2 = 1 / params_inc[["scale"]],
        width = width
        #M = length(si_vec),
        #y_vec = si_vec
      ),
      chains = 2,
      iter = 2000,
      seed = 42,
      verbose = TRUE
    )
    ## Save it now in case R crashes
    saveRDS(fit_4a_recall, glue::glue("stanfits/{prefix}{index}.rds"))
    fit_4a_recall
  }
)
