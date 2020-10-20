rstan::expose_stan_functions(here::here("stan-models/likelihoods.stan"))
set.seed(42)
nsim <- 100
alpha_invalid <- 1.5
beta_invalid <- 6.5
pinvalid <- 0.03

valid_si <- simulate_si(
  mean_inc, sd_inc, 6.5, 1.5, max_shed, nsim = nsim
)



invalid_si <- max(valid_si$si) *
  rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)

toss <- runif(nsim)

simulated_si <- c(
  invalid_si[toss <= pinvalid], valid_si$si[toss > pinvalid]
)

grid <- expand.grid(
  pinvalid = pinvalid,
  alpha1 = seq(0.1, 200, 1),
  beta1 = seq(0.1, 200, 1)
)

prior <- dunif(grid$pinvalid, log = TRUE) +
  dunif(grid$alpha1, alpha_invalid, 10, log = TRUE) +
  dunif(grid$beta1, 0, 10, log = TRUE)

alpha2 <- params_inc[["shape"]]
beta2 <- 1 / params_inc[["scale"]]
width <- min(simulated_si) / 2
max_si <- max(simulated_si) + 0.001

grid$posterior <- pmap_dbl(
  grid,
  function(pinvalid, alpha1, beta1) {
    ll <- sapply(
      simulated_si,
      function(x) {
        ll_valid <- scenario1a_lpdf(x, max_shed, alpha1, beta1,
                                    alpha2, beta2, width)
        ll_invalid <- dbeta(
          x/max_si, alpha_invalid, beta_invalid, log = TRUE
        )
        ll_valid + ll_invalid
      }
    )
    sum(ll)
  }
)

ggplot(grid, aes(alpha1, beta1, z = posterior)) +
  geom_contour_filled() +
  facet_wrap(~pinvalid)
