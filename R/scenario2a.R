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

simulated_2 <- simulated_2[simulated_2$t_1 < simulated_2$nu, ]

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
  chains = 3,
  iter = 10000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fitted_params <- rstan::extract(fits_2a)

x <- hermione::beta_shape1shape22muvar(
  fitted_params[["alpha1"]], fitted_params[["beta1"]]
)

x[["mu"]] <- max_shed * x[["mu"]]
x[["sigma2"]] <- max_shed^2 * x[["sigma2"]]
x[["sd"]] <- sqrt(x[["sigma2"]])

p1 <- ggplot(NULL, aes(x[["mu"]])) +
  geom_histogram(alpha = 0.3) +
  geom_vline(xintercept = mean_inf, linetype = "dashed")

p2 <- ggplot(NULL, aes(x[["sd"]])) +
  geom_histogram(alpha = 0.3) +
  geom_vline(xintercept = sd_inf, linetype = "dashed")

p <- cowplot::plot_grid(p1, p2, ncol = 1)

ggsave("infectious_profile_params.png", p)

nsamples <- length(fitted_params[[1]])
idx <- sample(nsamples, size = ceiling(nsamples/2), replace = FALSE)
shape1 <- fitted_params[[1]][idx]
shape2 <- fitted_params[[2]][idx]

si_post <- simulate_si(mean_inc, sd_inc, shape1, shape2, max_shed, 2, 2)

psi <- ggplot() +
  geom_density(
    data = simulated_2, aes(si),
    alpha = 0.3, fill = "blue", binwidth = 0.5
  ) +

  geom_density(
    data = si_post, aes(si),
    alpha = 0.3, fill = "red", binwidth = 0.5
  ) +
  geom_vline(
    xintercept = mean(simulated_2$si), col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Serial Interval")

ggsave("posterior_serial_interval.png", psi)
