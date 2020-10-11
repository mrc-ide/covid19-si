gen_data <- function(seed) {
  set.seed(seed + 16)
  list(
    N = 30, max_shed = 21,  width = 0.1
  )
}

gen_params <- function(seed, data) {
  set.seed(seed + 5e6)
  ##
  alpha1 <- runif(1, min = 1, max = 30)
  beta1 <- runif(1, min = 1, max = 30)
  message("alpha1 = ", alpha1)
  message("beta1 = ", beta1)
  alpha2 <- runif(1, min = 1, max = 30)
  beta2 <- runif(1, min = 1, max = 30)
  message("alpha2 = ", alpha2)
  message("beta2 = ", beta2)

  list(alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2)
}

gen_modeled_data <- function(seed, data, params) {
  set.seed(seed + 396)
  t_1 <- data$max_shed * stats::rbeta(
    n = data$N, shape1 = params$alpha1, shape2 = params$beta1
  )
  t_2 <- stats::rgamma(
    n = data$N, shape = params$alpha2, rate = params$beta2
  )
  si <- t_1 + t_2

  list(si = si)
}

sample_from_model <- function(seed, data, params, modeled_data, iters) {
  data$width <- min(modeled_data$si) / 2
  data_for_stan <- c(data, modeled_data)
  rstan::sampling(s1b_model, data = data_for_stan, seed = seed,
                  chains = 1, iter = 2 * iters, warmup = iters,
                  open_progress = FALSE, show_messages = FALSE,
                  refresh = 1000)
}

s1b_model <- stan_model(
  file = here::here("stan-models/scenario1b.stan")
)
set.seed(42)

s1b_sbc <- SBC$new(data = gen_data,
                   params = gen_params,
                   modeled_data = gen_modeled_data,
                   sampling = sample_from_model
                   )

##doParallel::registerDoParallel(cores = parallel::detectCores())
cal <- s1b_sbc$calibrate(N = 300, L = 15, keep_stan_fit = FALSE)
cal$summary()
cal$plot()


rank_a1 <- map_dbl(cal$calibrations, ~ .$ranks$alpha1)
hist(rank_a1, breaks = seq(0, max(rank_a1), by = 1))

rank_b1 <- map_dbl(cal$calibrations, ~ .$ranks$beta1)
hist(rank_b1, breaks = seq(0, max(rank_b1), by = 1))

rank_a2 <- map_dbl(cal$calibrations, ~ .$ranks$alpha2)
hist(rank_a2, breaks = seq(0, max(rank_a1), by = 1))

rank_b2 <- map_dbl(cal$calibrations, ~ .$ranks$beta2)
hist(rank_b2, breaks = seq(0, max(rank_b1), by = 1))
