f <- function(x, nu, alpha1, beta1, alpha2, beta2) {

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

g <- Vectorize(f)

simulated_2 <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  mean_inf = mean_inf,
  sd_inf = sd_inf,
  mean_iso = 2,
  sd_iso = 2,
  nsim = 1000
)

grid <- expand.grid(
  alpha1 = seq(1, 100, by = 1), beta1 = seq(0.1, 100, by = 0.5)
)

out <- pmap(
  grid,
  function(alpha1, beta1) {
    y <- g(simulated_2$si, simulated_2$nu, alpha1, beta1,
      params_inc[["shape"]], 1 / params_inc[["scale"]]
      )
    sum(y)
  }
)

simulated_2 <- simulated_2[simulated_2$t_1 < simulated_2$nu, ]

fits_2a <- stan(
  file = here::here("stan-models/scenario2a.stan"),
  data = list(
    N = nrow(simulated_2),
    si = simulated_2$si,
    nu = simulated_2$nu,
    max_shed = 21,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]]
  ),
  chains = 1,
  iter = 200,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


