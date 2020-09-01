simulated_1 <- simulate_si(
  mean_ip = mean_inc,
  sd_ip = sd_inc,
  mean_inf = mean_inf,
  sd_inf = sd_inf,
  nsim = 3000
)

nsim <- seq(from = 10, to = 50, by = 10)
names(nsim) <- nsim
ndatasets <- 10

sim_data <- map(
  nsim,
  function(n) {
    map(
      seq_len(ndatasets),
      function(i) {
        idx <- sample(seq_len(nrow(simulated_1)), n)
        out <- simulated_1[idx, ]
        out
      }
    )
  }
)


fits_1a <- map_depth(
  sim_data,
  2,
  function(df) {
      stan(
        file = here::here("stan-models/scenario1a.stan"),
        data = list(
          N = nrow(df),
          si = df$si,
          max_shed = 21,
          alpha2 = params_inc[["shape"]],
          beta2 = 1 / params_inc[["scale"]]
        ),
        chains = 3,
        iter = 6000
        ##control = list(adapt_delta = 0.99)
      )
  }
)
## Check convergence etc, using ggmcmc
## test_fit <- ggmcmc(ggs(fit_1a), here::here("figures/1a.pdf"))

## extract fits to turn alpha and beta into mu and cv
fitted_params <- map_depth(fits_1a, 2, rstan::extract)

fitted_params <- map_depth(
  fitted_params, 2, function(out) {
    x <- gamma_shapescale2mucv(
      out[["alpha1"]], 1 / out[["beta1"]]
    )
    out[["mu"]] <- x[["mu"]]
    out[["cv"]] <- x[["cv"]]
    out
  }
)

fitted_summary <- map_depth(
  fitted_params, 2, ~ map_dfr(., quantile, .id = "param")
)

fitted_summary <- map_dfr(
  fitted_summary, ~ bind_rows(., .id = "sim"), .id = "size"
)


### Figures
pmu <- ggplot() +
  geom_linerange(
    data = fitted_summary[fitted_summary$param %in% c("mu", "cv"), ],
    aes(x = sim, ymin = `25%`, ymax = `75%`),
  ) +
  facet_grid(size ~ param) +
  ##geom_vline(xintercept = mean_inf, col = "red", linetype = "dashed") +
  theme_minimal()
  ##xlab("Mean Infectious Period")


pcv <- ggplot() +
  geom_histogram(
    data = NULL, aes(x = params_inf_est[["cv"]], y = ..density..),
    alpha = 0.3
  ) +
  geom_vline(
    xintercept = sd_inf / mean_inf, col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("CV Infectious Period")

## Simulate with draws from posterior
posterior_mean <- sample(params_inf_est[["mu"]], size = 100)
posterior_sd <- sample(
  params_inf_est[["mu"]] * params_inf_est[["cv"]], size = 100
)

simulated_post <- map2_dfr(
  posterior_mean,
  posterior_sd,
  function(x, y) simulate_si(mean_inc, sd_inc, x, y, nsim = 20)
)

psi <- ggplot() +
  geom_histogram(
    data = simulated_1, aes(x = si, y = ..density..),
    alpha = 0.3, fill = "blue"
  ) +

  geom_histogram(
    data = simulated_post, aes(x = si, y = ..density..),
    alpha = 0.3, fill = "red"
  ) +
  geom_vline(
    xintercept = mean(simulated_1$si), col = "red", linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Serial Interval")









