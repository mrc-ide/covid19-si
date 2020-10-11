gen_data <- function(seed) {
  set.seed(seed + 1e6)
  list(
    N = 30, max_shed = 21,
    alpha2 = params_inc$shape,
    beta2 = 1 / params_inc$scale
  )
}

gen_params <- function(seed, data) {
  set.seed(seed + 2e6)
  alpha1 <- runif(1, min = 1, max = 50)
  beta1 <- runif(1, min = 1, max = 50)
  list(alpha1 = alpha1, beta1 = beta1)
}

gen_modeled_data <- function(seed, data, params) {
  set.seed(seed + 3e6)
  t_1 <- data$max_shed * stats::rbeta(
    n = data$N, shape1 = params$alpha1, shape2 = params$beta1
  )
  t_2 <- stats::rgamma(
    n = data$N, shape = data$alpha2, rate = data$beta2
  )
  list(si = t_1 + t_2)
}

sample_from_model <- function(seed, data, params, modeled_data, iters) {
  data_for_stan <- c(data, modeled_data)
  rstan::sampling(s1a_model, data = data_for_stan, seed = seed,
                  chains = 1, iter = 2 * iters, warmup = iters,
                  open_progress = FALSE, show_messages = FALSE,
                  refresh = 1000)
}

s1a_model <- stan_model(
  file = here::here("stan-models/scenario1a.stan")
)

s1a_sbc <- SBC$new(data = gen_data,
                   params = gen_params,
                   modeled_data = gen_modeled_data,
                   sampling = sample_from_model
                   )

doParallel::registerDoParallel(cores = parallel::detectCores())
cal <- s1a_sbc$calibrate(N = 600, L = 30, keep_stan_fit = FALSE)

saveRDS(cal, "cal_1a.rds")


rank_a1 <- map_dbl(cal$calibrations, ~ .$ranks$alpha1)

p <- ggplot() +
  geom_histogram(data = NULL, aes(rank_a1), binwidth = 1, alpha = 0.2) +
  geom_hline(
    yintercept = qbinom(p = 0.025, 600, 1/31), linetype = "dashed"
  ) +
  geom_hline(
    yintercept = qbinom(p = 0.975, 600, 1/31), linetype = "dashed"
  ) +
  theme_minimal() + xlab("shape1")


ggsave("figures/hist_rank_a1_1a.png")

rank_b1 <- map_dbl(cal$calibrations, ~ .$ranks$beta)

p <- ggplot() +
  geom_histogram(data = NULL, aes(rank_b1), binwidth = 1, alpha = 0.2) +
  geom_hline(
    yintercept = qbinom(p = 0.025, 600, 1/31), linetype = "dashed"
  ) +
  geom_hline(
    yintercept = qbinom(p = 0.975, 600, 1/31), linetype = "dashed"
  ) +
  theme_minimal() + xlab("shape2")

ggsave("figures/hist_rank_b1_1a.png")
