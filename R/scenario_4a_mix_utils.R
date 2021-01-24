simulate_4a_mix <- function(params_inc, params_inf, params_iso,
                            offset, max_shed, pinvalid, nsim,
                            alpha_invalid, beta_invalid, min_si, max_si) {

  ## Simulate lots more than you need because a lot will get filtered
  ## out when you do t_1 < nu
  valid_si <- better_simulate_si(
    params_inc, params_inf, params_iso,
    min_inf = offset, max_inf = max_shed, nsim * 10
  )
  unconditional <- valid_si
  valid_si <- valid_si[valid_si$t_1 < valid_si$nu, ]
  ## Make sure you have enough else there will be NAs
  message("Need ", nsim, " valid observations. Have ", nrow(valid_si))
  idx <- sample(nrow(valid_si), nsim, replace = TRUE)
  valid_si <- valid_si[idx, ]

  invalid_si <- rbeta(
    nsim, shape1 = alpha_invalid, shape2 = beta_invalid
  )
  f <- map_into_interval(0, 1, min_si, max_si)
  invalid_si <- f(invalid_si)
  invalid_nu <- stats::rgamma(
    n = nsim, shape = params_iso$shape, scale = params_iso$scale
  )
  invalid_si <- data.frame(si = invalid_si, nu = invalid_nu)

  toss <- runif(nsim)
  sim_valid <- data.frame(
    valid_si[toss > pinvalid, c("si", "nu")], group = "valid"
  )
  if (! any(toss <= pinvalid)) {
    ## This will happen if pinvalid = 0
    ##sim_invalid <- data.frame(si = NA, nu = NA, group = "invalid")
    sim_invalid <- NULL
    simulated_si <- sim_valid
  } else {
    sim_invalid <- data.frame(
      si = invalid_si[toss <= pinvalid, ], group = "invalid"
    )
    simulated_si <- rbind(sim_invalid, sim_valid)
  }

  list(
    valid_si = valid_si, invalid_si = invalid_si,
    simulated_si = simulated_si, unconditional = unconditional
  )
}
