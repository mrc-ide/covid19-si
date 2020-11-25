## Returns the pdf, you then intergrate it to get the constant
full_model <- function(x, nu, max_shed, offset1, recall,
                       alpha1, beta1, alpha2, beta2, width) {
    if (x > max_shed) ulim <- max_shed
    else ulim <- x
    if(ulim > nu) ulim <- nu
    out <- 0
    max_shed_shifted <- max_shed - offset1
    nu_shifted <- nu - offset1
    s <- offset1 + width
    while (s < ulim) {
      inf_density <- dbeta(
      (s - offset1)/max_shed_shifted, alpha1, beta1, log = TRUE
      )
      inc_density <- dgamma(x - s, alpha2, beta2, log = TRUE)
      out <- out + exp(inf_density + inc_density)
      s = s + width
    }
    out <- log(out)
    if(nu < max_shed) {
      out <- out -
        pbeta(nu_shifted / max_shed_shifted, alpha1, beta1, log = TRUE)
    }
    out <- out - recall * abs(x - nu)
    out
}


quantile_as_df <- function(vec, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {

  out <- data.frame(val = quantile(vec, probs), check.names = FALSE)
  out <- tibble::rownames_to_column(out, "var")
  out <- rbind(
    out,
    data.frame(var = "mu", val = mean(vec)),
    data.frame(var = "sd", val = sd(vec))
  )
  out
}

## More aligned args, and also allows for pre-symptpmatic infection
better_simulate_si <- function(params_inc, params_inf, params_iso,
                               min_inf, max_inf, nsim = 50) {

  t_1 <- stats::rbeta(
    n = nsim, shape1 = params_inf$shape1, shape2 = params_inf$shape2
    )
  ## Map the possible infection times into interval
  ## (min_inf, max_inf)
  f <- map_into_interval(0, 1, min_inf, max_inf)
  t_1 <- f(t_1)

  t_2 <- stats::rgamma(
    n = nsim, shape = params_inc$shape, rate = 1 / params_inc$scale
  )

  out <- data.frame(t_1 = t_1, t_2 = t_2, si = t_1 + t_2)


  out$nu <- stats::rgamma(
      n = nsim, shape = params_iso$shape, rate = 1 / params_iso$scale
    )

  out
}

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

# beta mu and cv <-> shape 1 and 2


beta_muvar2shape1shape2 <- function(mu, sigma2) {

  shape1 <- (mu^2 * (1 - mu) / sigma2) - mu
  shape2 <- shape1 * (1 - mu) / mu
  list(shape1 = shape1, shape2 = shape2)
}

beta_shape1shape22muvar <- function(shape1, shape2) {

  mu <- shape1 / (shape1 + shape2)
  sigma2 <- (shape1 * shape2) / ((shape1 + shape2)^2 * (shape1 + shape2 + 1))
  list(mu = mu, sigma2 = sigma2)
}

## Linear map of one interval into another
## maps (x_1, y_1) into (x_2, y_2)
## Returns a function
map_into_interval <- function(x_1, y_1, x_2, y_2) {

  slope <- (x_2 - y_2) / (x_1 - y_1)
  intercept <- x_2 - x_1 * slope
  f <- function(x) slope * x + intercept
  f
}
