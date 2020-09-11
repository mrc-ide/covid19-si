simulate_si <- function(mean_ip,
                        sd_ip,
                        shape1_inf,
                        shape2_inf,
                        max_shed,
                        mean_iso = NULL,
                        sd_iso = NULL,
                        nsim = 50) {
  params_ip <- epitrix::gamma_mucv2shapescale(
    mu = mean_ip, cv = sd_ip / mean_ip
  )

  t_1 <- max_shed * stats::rbeta(
    n = nsim, shape1 = shape1_inf, shape2 = shape2_inf
  )

  t_2 <- stats::rgamma(
    n = nsim, shape = params_ip$shape, rate = 1 / params_ip$scale
  )

  out <- data.frame(t_1 = t_1, t_2 = t_2, si = t_1 + t_2)

  if (!is.null(mean_iso)) {

    params_iso <- epitrix::gamma_mucv2shapescale(
      mu = mean_iso, cv = sd_iso / mean_iso
    )
    out$nu <- stats::rgamma(
      n = nsim, shape = params_iso$shape, rate = 1 / params_iso$scale
    )
  }
  ## 0 Serial Interval is not allowed if we do not account for
  ## pre-symptomatic infectivity.


  out
}

## Scenario 2 log-likelihood as implemented in Stan.
scenario2_ll <- function(x, nu, alpha1, beta1, alpha2, beta2) {

  inf_density <- dgamma(
    0.5, shape = alpha1, rate = beta1, log = TRUE
  ) + log(0.5)
  inc_density <- dgamma(
    0.5, shape = alpha2, rate = beta2, log = TRUE
  ) + log(0.5)
  out <- exp(inf_density + inc_density)
  s <- 0.5
  while (s <= x) {
    s <- s + 0.5
    inf_density <- dgamma(
      0.5, shape = alpha1, rate = beta1, log = TRUE
    ) + log(0.5)
    inc_density <- dgamma(
      0.5, shape = alpha2, rate = beta2, log = TRUE
    ) + log(0.5)
    out <- out + exp(inf_density + inc_density)
  }
  out <- log(out)
  out <- out - pgamma(nu, alpha1, beta1, log = TRUE)
  out
}

vscenario2_ll <- Vectorize(scenario2_ll)
