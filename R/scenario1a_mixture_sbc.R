gen_data <- function(seed) {
  set.seed(seed + 16)
  list(
    N = 30, max_shed = 21, alpha2 = params_inc$shape,
    beta2 = 1 / params_inc$scale, alpha_invalid = 0.5,
    beta_invalid = 0.5, max_si = 40, min_si = -2, width = 0.1
  )
}

gen_params <- function(seed, data) {
  set.seed(seed + 5e6)
  ##
  pinvalid <- runif(1, min = 0, max = 1)
  ##theta <- c(pinvalid, 1 - pinvalid)
  alpha1 <- runif(1, min = data$alpha_invalid, max = 100)
  beta1 <- runif(1, min = 0, max = 100)
  message("pinvalid = ", pinvalid)
  message("alpha1 = ", alpha1)
  message("beta1 = ", beta1)
  list(pinvalid = pinvalid, alpha1 = alpha1, beta1 = beta1)
}

gen_modeled_data <- function(seed, data, params) {
  set.seed(seed + 396)
  t_1 <- data$max_shed * stats::rbeta(
    n = data$N, shape1 = params$alpha1, shape2 = params$beta1
  )
  t_2 <- stats::rgamma(
    n = data$N, shape = data$alpha2, rate = data$beta2
  )
  valid_si <- t_1 + t_2

  invalid_si <- (data$max_si - data$min_si) *
    rbeta(data$N, data$alpha_invalid, data$beta_invalid)
  toss <- runif(data$N)
  si <- c(
    invalid_si[toss <= params$pinvalid],
    valid_si[toss > params$pinvalid]
  )
  list(si = si)
}

sample_from_model <- function(seed, data, params, modeled_data, iters) {
  data$width <- min(modeled_data$si[modeled_data$si > 0]) / 2
  data$max_si <- max(modeled_data$si) + 0.001
  data$min_si <- min(modeled_data$si) - 0.001
  data_for_stan <- c(data, modeled_data)
  iters <- 2000
  rstan::sampling(mixture_model, data = data_for_stan, seed = seed,
                  chains = 2, iter = 2 * iters, warmup = iters,
                  open_progress = FALSE, show_messages = FALSE,
                  refresh = 1000)
}

mixture_model <- stan_model(
  file = here::here("stan-models/scenario1a_mixture_general.stan")
)
set.seed(42)
mixture_sbc <- SBC$new(data = gen_data,
                       params = gen_params,
                       modeled_data = gen_modeled_data,
                       sampling = sample_from_model
                       )

##doParallel::registerDoParallel(cores = parallel::detectCores())
cal <- mixture_sbc$calibrate(N = 10, L = 15, keep_stan_fit = TRUE)
cal$summary()
cal$plot()

cal$plot("beta1")

rank_p <- map_dbl(cal$calibrations, ~ .$ranks$pinvalid)
hist(rank_p, breaks = seq(0, max(rank_p), by = 1))

rank_a1 <- map_dbl(cal$calibrations, ~ .$ranks$alpha1)
hist(rank_a1, breaks = seq(0, max(rank_a1), by = 1))

rank_b1 <- map_dbl(cal$calibrations, ~ .$ranks$beta)
hist(rank_b1, breaks = seq(0, max(rank_b1), by = 1))

#######################
######## Could the approximation be an issue?
expose_stan_functions(stanmodel = mixture_model)

ll_summed1 <- map_dfr(
  seq(0.001, 2, by = 0.001),
  function(w) {
    map_dfr(
      1:100,
      function(x) {
        out <- scenario1a_lpdf(
          x, 21, alpha1 = params_inf$shape1, beta1 = params_inf$shape2,
          alpha2 = params_inc$shape, beta2 = 1 / params_inc$scale,
          width = w
        )
        data.frame(x = x, ll = out, width = w)
      }
    )
  }
)

ll_integrated <- map_dfr(
  seq(0.001, 2, by = 0.001),
  function(w) {
    map_dfr(
      1:100,
      function(x) {
        f <- function(s) {
          dbeta(
            s / 21, shape1 = params_inf$shape1, shape2 = params_inf$shape2
          ) * dgamma(
                x - s, shape = params_inc$shape, rate = 1 / params_inc$scale
              )
        }
        out <- stats::integrate(f, 0, x, stop.on.error = FALSE)
        data.frame(x = x, ll = log(out$value), width = w)
      }
    )
  }
)

error <- data.frame(
  x = ll_summed1$x,
  width = ll_summed1$width,
  err = ll_summed1$ll - ll_integrated$ll
)


ggplot(error) +
  geom_point(aes(width,  err, col = x))
  scale_fill_distiller(palette = "Greens")
