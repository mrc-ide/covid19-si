## Set up params grid
param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2", "inf_par3"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = c("iso_par1", "iso_par2", "iso_par3"),
  params_offset = c("offset1", "offset2", "offset3"),
  params_recall = c("recall1", "recall2", "recall3"),
  params_pinv = c("pinvalid1", "pinvalid2", "pinvalid3"),
  stringsAsFactors = FALSE
)

param_grid <- param_grid[1, ]
nsim <- 20000

params_inf_all <- map(
  param_grid$params_inf,
  function(x) {
    out <- params[[x]]
    beta_muvar2shape1shape2(
      out$mean_inf/max_shed, out$sd_inf^2 / max_shed^2
    )
  }
)

params_inc_all <- map(
  param_grid$params_inc,
  function(x) {
    out <- params[[x]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_iso_all <- map(
  param_grid$params_iso,
  function(x) {
    out <- params[[x]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_pinv_all <- 0

params_offset_all <- map(
  param_grid$params_offset, function(x) params[[x]]
)

params_recall_all <- map(
  param_grid$params_recall, function(x) params[[x]]
)

prefix <- "full_model_sim_"
nsim_post_filter <- 100

simulated_data <- pmap(
  list(
    params_inf = params_inf_all,
    params_inc = params_inc_all,
    params_iso = params_iso_all,
    params_offset = params_offset_all
  ),
  function(params_inf, params_inc, params_iso, params_offset) {
    sim_data <- better_simulate_si(
      params_inc, params_inf, params_iso, params_offset, max_shed, nsim
    )
    sim_data <- sim_data[sim_data$t_1 <= sim_data$nu, ]
    sim_data <- sim_data[abs(sim_data$si) > 0.1, ]
    ## Make sure we have at least 200 rows.
    idx <- sample(nrow(sim_data), nsim_post_filter, replace = TRUE)
    sim_data[idx, ]
  }
)

invalid_si <- map(
  params_iso_all,
  function(params_iso) {
    invalid_si <- rbeta(
      nsim_post_filter, shape1 = alpha_invalid, shape2 = beta_invalid
    )

    invalid_iso <- rgamma(
      nsim_post_filter, shape = params_iso$shape, scale = params_iso$scale
    )
    data.frame(si = invalid_si, nu = invalid_iso)
  }
)


invalid_remapped <- map2(
  simulated_data,
  invalid_si,
  function(valid, invalid) {
    max_si <- max(valid$si)
    ## invalid SIs are draws from beta. Map them into
    ## min and max of valid SI
    ##
    f <- map_into_interval(0, 1, 0.5 * min(valid$si), 2 * max(valid$si))
    invalid$si <- f(invalid$si)
    invalid
  }
)

mixed <- pmap(
  list(
    valid = simulated_data,
    invalid = invalid_remapped,
    pinvalid = params_pinv_all
  ),
  function(valid, invalid, pinvalid) {
    ##pinvalid <- params[[params_pinv]]
    toss <- runif(nrow(valid))
    valid$type <- "valid"
    invalid$type <- "invalid"
    rbind(
      valid[toss > pinvalid , c("si", "nu", "type")],
      invalid[toss <= pinvalid ,c("si", "nu", "type")]
    )
  }
)

denominator <- function(nu, recall, max_si, min_si) {
  if (nu > min_si & nu < max_si) {
    out <- (2 - exp(recall * (-nu + min_si)) -
              exp(recall * (-max_si + nu)))
  } else if (nu >= max_si) {
    out <- exp(recall * (max_si - nu)) - exp(recall * (min_si - nu))
  } else if (nu <= min_si) {
    out <- -exp(recall * (-max_si + nu)) + exp(recall * (-min_si + nu))
  }
  out / recall
}

f <- Vectorize(denominator)

with_recall_bias <- pmap(
  list(
    df = mixed,
    recall_true = params_recall_all,
    offset = params_offset_all
  ),
  function(df, recall_true, offset) {

    ##x <- f(df$nu, recall_true, max(df$si), offset)
    df$p_si <- exp(abs(df$si - df$nu) * -recall_true)
    idx <- sample(nrow(df), nrow(df), replace = TRUE, prob = df$p_si)
    df[idx, ]
  }
)


width <- 0.1

fits <- pmap(
  list(
    sim_data = with_recall_bias,
    param_inc = params_inc_all,
    param_offset = params_offset_all
  ),
  function(sim_data, param_inc, param_offset) {

    fit <- stan(
      file = here::here("stan-models/full_model.stan"),
      data = list(
        N = length(sim_data$si),
        si = sim_data$si,
        nu = sim_data$nu,
        max_shed = max_shed,
        offset1 = param_offset,
        max_si = max(sim_data$si) + 0.001,
        min_si = param_offset, ## assuming the smallest incubation period is 0
        alpha2 = param_inc[["shape"]],
        beta2 = 1 / param_inc[["scale"]],
        ##alpha_invalid = alpha_invalid,
        ##beta_invalid = beta_invalid,
        width = width
      ),
      chains = 3,
      iter = 2000,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
  }
)


grid <- expand.grid(recall = seq(0, 2, by = 0.01))

ll <- pmap_dfr(
  grid,
  function(recall) {
    pmap_dfr(
      sim_data[, c("si", "nu")],
      function(si, nu) {
        alpha1 <- params_inf_all[[1]][[1]]
        beta1 <- params_inf_all[[1]][[2]]
        alpha2 <- params_inc_all[[1]][[1]]
        beta2 <- params_inc_all[[1]][[2]]
        out <- full_model_lpdf(
          si, nu, max_shed, -1, recall,
          alpha1, beta1, alpha2,
          1 / beta2, 0.01,
          max(sim_data$si) + 0.001, -1
        )
        data.frame(
          si = si, nu = nu, ll = out, recall = recall
        )
      }
    )
  }
)

x <- group_by(ll, recall) %>%
  summarise(total = sum(ll)) %>%
  ungroup()

x[which.max(x$total), ]

total_ll <- map_dbl(ll, sum)
ggplot() + geom_point(aes(grid$recall, total_ll))
## grid <- expand.grid(si = -1:40, nu = 0:40)
## den <-  map_dbl(grid$nu, ~ denominator(., 0.1, 40, -1))
## num <- exp(-0.1 * abs(grid$nu - grid$si))
## grid$p_si <- num / den


out1 <- pmap(
  sim_data[, c("si", "nu")],
  function(si, nu) {
    full_model_lpdf(si, nu, max_shed, -1, 0.01,
                    3.552653, 6.896327, param_inc[["shape"]],
                    1 / param_inc[["scale"]], width,
                    max(sim_data$si) + 0.001, -1)
    ## invalid_lpdf(si,  max(sim_data$si) + 0.001, param_offset,
    ##              alpha_invalid, beta_invalid)
  }
)


out2 <- pmap(
  sim_data[, c("si", "nu")],
  function(si, nu) {
    full_model_lpdf(si, nu, max_shed, -1, 0.01,
                    7, 10, param_inc[["shape"]],
                    1 / param_inc[["scale"]], width,
                    max(sim_data$si) + 0.001, -1)
  }
)

out3 <- pmap(
  sim_data[, c("si", "nu")],
  function(si, nu) {
    full_model_lpdf(si, nu, max_shed, -1, 0.01,
                    70, 10, param_inc[["shape"]],
                    1 / param_inc[["scale"]], width,
                    max(sim_data$si) + 0.001, -1)
  }
)


ggplot() +
  geom_line(aes(sim_data$si, unlist(out1)), col = "red") +
  geom_line(aes(sim_data$si, unlist(out2)), col = "blue") +
  geom_line(aes(sim_data$si, unlist(out3)))



pmap(
  sim_data[1:2, c("si", "nu")],
  function(si, nu) {
    alpha1 <- params_inf_all[[1]][[1]]
    beta1 <- params_inf_all[[1]][[2]]
    alpha2 <- params_inc_all[[1]][[1]]
    beta2 <- params_inc_all[[1]][[2]]
    full_model_lpdf(x = si, nu = nu, max_shed = max_shed, offset1 = -1,
                    width = 0.01, alpha1, beta1,
                    alpha2, 1/ beta2,
                    width, 0.001 + max(sim_data$si), -1
                    )
  })

normalising_constant(real nu, real max_si, real min_si,
                     real recall)


denominator(nu = 3.87151, recall = 0.01, max_si = 21.1471, min_si = -1)
nu <- seq(0, 5, by = 0.1)
alpha1 <- params_inf_all[[1]][[1]]
beta1 <- params_inf_all[[1]][[2]]
alpha2 <- params_inc_all[[1]][[1]]
beta2 <- params_inc_all[[1]][[2]]

out <- map_dbl(
  nu, function(x) {
    normalising_constant(x, max_shed, -1, 0.5,
                         alpha1, beta1, alpha2, beta2, width, 40, -1)
  }
)


out2 <- map_dbl(
  sim_data$nu, ~ log(denominator(
                nu = ., max_si = max(sim_data$si) + 0.001,
                min_si = -1, recall = 0.01))
)
