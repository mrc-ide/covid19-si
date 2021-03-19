library(context)
root <- "20210319"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue")
source_files <- c("global.R", "utils.R")
ctx <-context_save(
  root, packages = packages, sources = source_files
)

###
prefix <- "4a_mix_with_normalisation_sim"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = "iso_par1",
  params_offset = c("offset1", "offset2"),
  params_pinvalid = c("pinvalid1", "pinvalid2"),
  stringsAsFactors = FALSE
)


index <- 1
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

## with -ve nu
unconditional_data <- pmap(
  list(
   dat = uncdtnl_data,
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

## conditional
simulated_data <- map(
  unconditional_data,
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

outfiles <- glue::glue("data/{prefix}_{seq_along(sampled)}_data.rds")
walk2(sampled, outfiles, function(x, y) saveRDS(x, y))

obj <- didehpc::queue_didehpc(ctx)

grp <- obj$lapply(
             index, function(i) {
               params_inc <- params_inc_all[[i]]
               params_offset <- params_offsets_all[[i]]
               params_iso <- params_iso_all[[i]]
               sim_data <- sampled[[i]]
               si_vec <- seq(params_offset + 0.5, max_valid_si, 1)
               width <- 0.1
               sim_data <- dplyr::arrange(sim_data, nu)
               fit_4a <- stan(
                 file = "scenario4a_mixture.stan",
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
                 seed = 42,
                 verbose = TRUE
               )
               outfile <- glue::glue("stanfits/{prefix}_{index}.rds")
               saveRDS(fit_4a, outfile)
               fit_4a
             }
           )
