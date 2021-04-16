source("global2.R")
source("utils.R")
prefix <- "full_model"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = c("iso_par1", "iso_par2"),
  params_offset = c("offset1", "offset2", "offset3"),
  params_pinvalid = c("pinvalid1", "pinvalid2", "pinvalid3"),
  params_beta = c("recall1", "recall2", "recall3"),
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
params_offsets_all <- map(
  param_grid$params_offset,
  function(params_offset) params[[params_offset]]
)

params_pinv <- map(
  param_grid$params_pinv, function(params_pinv) params[[params_pinv]]
)

params_beta <- map(
  param_grid$params_beta, function(params_beta) params[[params_beta]]
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

outfiles <- glue::glue("data/{prefix}_{seq_along(mixed)}_data.rds")
walk2(sampled, outfiles, function(x, y) saveRDS(x, y))
