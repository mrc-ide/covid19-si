vscenario_1 <- Vectorize(probability_basic)

ll_scenario1 <- function(t, inf_distr, inc_distr, inf_params, inc_params) {
  message(paste(inf_params, collapse = " "))
  out <- vscenario_1(t, inf_distr, inc_distr, inf_params, inc_params)
  sum(log(out))
}

mh_scenario1 <- function(t, start, n_iter, prop_sd, params_prior, inc_params) {

  chain <- matrix(NA, ncol = length(start), nrow = n_iter + 1)
  posterior <- matrix(NA, ncol = 1, nrow = n_iter + 1)
  chain[1, ] <- start
  ## Uniform prior
  prior <- sum(dunif(t, min(t), max(t), log = TRUE))
  inf_params <- list(shape = start[1], scale = start[2])
  ll <- ll_scenario1(t, "dgamma", "dgamma", inf_params, inc_params)
  posterior[1] <- ll + prior

  for (idx in 1:n_iter){

    proposal <- chain[idx, ] + rnorm(length(start), 0, prop_sd)
    if (any(proposal < 0)) {
      chain[idx + 1,] <- chain[idx, ]
      posterior[idx + 1] <- posterior[idx]
    } else {
      inf_params <- list(shape = proposal[1], scale = proposal[2])
      posterior_new <- prior +
        ll_scenario1(
          t, "dgamma", "dgamma", list(shape = proposal[1], scale = proposal[2]), inc_params
      )
      posterior_old <- prior +
        ll_scenario1(
          t, "dgamma", "dgamma", list(shape = chain[idx, 1], scale = chain[idx, 2]), inc_params
        )
      prob <- exp(posterior_new - posterior_old)
      if (runif(1) < prob) {
        chain[idx + 1, ] <- proposal
        posterior[idx + 1] <- posterior
      }
      else {
        chain[idx + 1,] <- chain[idx, ]
        posterior[idx + 1] <- posterior[idx]
      }
    }
  }
  list(
    chain = chain, posterior = posterior
  )
}




inf_params <- list(mu = 2.2, cv = 3 / 2.2)
inc_params <- list(mu = 6.48, cv = 3.83 / 6.48)

x <- hermione::simulate_si(
  mean_ip = inc_params$mu,
  sd_ip = inc_params$mu * inc_params$cv,
  mean_inf = inf_params$mu,
  sd_inf = inf_params$mu * inf_params$cv,
  nsim = 20
)
out <- mh_scenario1(x$si, start, 10000, 1, NULL, inc_params)
