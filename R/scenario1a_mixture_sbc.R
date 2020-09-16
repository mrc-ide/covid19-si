gen_data <- function(seed) {
  set.seed(seed + 1e6)
  list(
    N = 100, max_shed = 21, alpha2 = 1, beta2 = 1, alpha_invalid = 1,
    beta_invalid = 1, max_si = 40, min_si = 0
  )
}

gen_params <- function(seed, data) {
  set.seed(seed + 2e6)
  pinvalid <- runif(1, min = 0, max = 1)
  alpha1 <- runif(1, min = 1, max = 100)
  beta1 <- runif(1, min = 1, max = 100)
  list(pinvalid = pinvalid, alpha1 = alpha1, beta1 = beta1)
}

gen_modeled_data <- function(seed, data, params) {
  set.seed(seed + 3e6)
  t_1 <- data$max_shed * stats::rbeta(
    n = data$N, shape1 = params$alpha1, shape2 = params$beta1
    )
  t_2 <- stats::rgamma(
    n = data$N, shape = data$alpha2, rate = data$beta2
  )
  valid_si <- t_1 + t_2

  invalid_si <- data$max_si *
    rbeta(data$N, data$alpha_invalid, data$beta_invalid)
  toss <- runif(data$N)
  si <- c(
    invalid_si[toss < params$pinvalid],
    valid_si[toss >= params$pinvalid]
  )
  list(si = si)
}

sample_from_model <- function(seed, data, params, modeled_data, iters) {
  data_for_stan <- c(data, modeled_data)
  rstan::sampling(mixture_model, data = data_for_stan, seed = seed,
                  chains = 1, iter = 2 * iters, warmup = iters,
                  open_progress = FALSE, show_messages = FALSE,
                  refresh = 1000)
}

mixture_model <- stan_model(
  file = here::here("stan-models/scenario1a_mixture.stan")
)

mixture_sbc <- SBC$new(data = gen_data,
                       params = gen_params,
                       modeled_data = gen_modeled_data,
                       sampling = sample_from_model)

doParallel::registerDoParallel(cores = parallel::detectCores())
mixture_sbc$calibrate(N = 32, L = 20, keep_stan_fit = FALSE)
