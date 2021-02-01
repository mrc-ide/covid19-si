## Set up params grid
prefix <- "4a_mix_with_normalisation_sim_"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = c("iso_par1", "iso_par2"),
  params_offset = c("offset1", "offset2", "offset3"),
  params_pinvalid = c("pinvalid1", "pinvalid2", "pinvalid3"),
  stringsAsFactors = FALSE
)

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
    simulate_4a_mix(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_pre_filter, alpha_invalid, beta_invalid,
      min_si, max_si
    )
  }
)

iwalk(
  simulated,
  function(df, index) {
    x <- df[["valid_si"]]
    y <- df[["unconditional"]]
    mixed <- df[["simulated_si"]]
    p <- ggplot() +
      geom_density(aes(x$si, fill = "red"), col = NA, alpha = 0.3) +
      geom_density(aes(y$si, fill = "blue"), col = NA, alpha = 0.3) +
      geom_density(aes(mixed$si, fill = "green"), col = NA, alpha = 0.3) +
      scale_fill_identity(
        breaks = c("red", "blue", "green"),
        labels = c("Conditional on nu", "Unconditional", "Mixed"),
        guide = "legend"
      ) +
      theme(legend.position = "top", legend.title = element_blank())
    ggsave(glue::glue("figures/{prefix}{index}_simulated.png"))
  }
)

sampled <- map(
  simulated, function(df) {
    ## might have NA from when pinvalid = 0
    out <- df[["simulated_si"]]
    out <- na.omit(out)
    ## Round here because if nu becomes 0 after rounding, we want to
    ## throw it out
    ##out$si <- round(out$si)
    ##out$nu <- round(out$nu)
    ##out <- out[out$nu > 0, ]

    idx <- sample(nrow(out), nsim_post_filter, replace = TRUE)
    out[idx, ]
  }
)


fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offset_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, sim_data, index) {
    ## Rounding now to check things
    min_si <- params_offset - 0.001
    max_si <- max(sim_data$si) + 1
    si_vec <- seq(min_si + 0.001, max_si, 1)
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
        max_si = max_si,
        min_si = min_si,
        width = width,
        M = length(si_vec),
        y_vec = si_vec
      ),
      ##chain = 1, iter = 1000,
      seed = 42,
      verbose = TRUE
    )
    ## Save it now in case R crashes
    saveRDS(fit_4a, glue::glue("stanfits/{prefix}{index}.rds"))
    fit_4a
  }
)

## infiles <- glue::glue("stanfits/{prefix}{index}.rds")
## fits <- map(infiles, readRDS)
## Debugging
## sim_data <- sampled[[1]]
## min_si <- 0.001
## si_vec <- seq(min_si, max_si, 1)
## out <- map_dbl(sim_data$nu, function(nu) {
##  s4_normalising_constant(si_vec, nu, 21, 0, 3.5, 1 / 33.5, 4, 2, 0, 40, 0.1)
## })
## rstan::expose_stan_functions("stan-models/likelihoods.stan")
## out <- pmap_dbl(sim_data[1:2, c("si", "nu")], function(si, nu) {
##  scenario4a_lpdf(si, nu, 21, 0, 3.5, 1 / 33.5, 4, 2, 0, 40, 0.1)
## })

## grid <- expand.grid(alpha1 = seq(1, 5, 0.5), beta1 = seq(1, 40, 0.5))

## valid <- pmap_dbl(
##   grid, function(alpha1, beta1) {
##     out <- pmap_dbl(
##       sim_data[, c("si", "nu")],
##       function(si, nu) {
##         scenario4a_lpdf(
##           si, nu, 21, 0, alpha1, beta1, 4, 2, 0.01, 10, 0.1
##         ) -  s4_normalising_constant(si_vec, nu, 21, 0, alpha1, beta1,
##                                      4, 2, 0.01, 10, 0.1)

##       }
##     )
##     sum(out)
##   }
## )


## invalid <- pmap_dbl(
##   grid, function(alpha1, beta1) {
##     out <- pmap_dbl(
##       sim_data[, c("si", "nu")],
##       function(si, nu) {
##         invalid_lpdf(si, min_si, max_si, alpha_invalid, beta_invalid)
##       }
##     )
##     sum(out)
##   }
## )

## grid$ll <- ll

## ggplot(grid, aes(alpha1, beta1, fill = ll)) +
##   geom_tile() +
##   scale_fill_distiller(palette = "Greens", direction = -1)
