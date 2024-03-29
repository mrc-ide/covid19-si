## params is a named list of estimated parameters
log_prior_nf <- function(params, mixture, recall, a_priors, b_priors) {
  prob_a <- dnorm(params$a, mean = a_priors$mean, sd = a_priors$sd) / pnorm(0, mean = a_priors$mean, sd = a_priors$sd)
  prob_b <-  dnorm(params$b, mean = b_priors$mean, sd = b_priors$sd) / pnorm(0, mean = b_priors$mean, sd = b_priors$sd)
  prob_c <- dunif(params$c, min = 0, max = 1)
  prob_tmax <- dunif(params$tmax, min = -20, max = 10)
  prob_pinvalid <- ifelse(mixture, dbeta(params$pinvalid, shape1 = 4, shape2 = 10), 1)
  prob_recall <- ifelse(recall, dunif(params$recall, min = 0, max = 5), 1)

  log(prob_a) + log(prob_b) + log(prob_c) + log(prob_tmax) +
    log(prob_pinvalid) + log(prob_recall)
}


log_prior_beta <- function(params, mixture, recall) {

  prob_a <- dunif(params$alpha1, min = 0, max = 100)
  prob_b <-  dunif(params$beta1, min = 0, max = 100)
  prob_pinvalid <- ifelse(mixture, dbeta(params$pinvalid, shape1 = 4, shape2 = 10), 1)
  prob_recall <- ifelse(recall, dunif(params$recall, min = 0, max = 5), 1)

  log(prob_a) + log(prob_b) +
    log(prob_pinvalid) + log(prob_recall)
}

log_prior_skew_normal <- function(params, mixture, recall) {

  prob_a <- dunif(params$a, min = -10, max = 10)
  prob_b <-  dunif(params$b, min = 0, max = 5)
  prob_c <-  dunif(params$c, min = -10, max = 10)
  prob_pinvalid <- ifelse(
    mixture, dbeta(params$pinvalid, shape1 = 4, shape2 = 10), 1
  )
  prob_recall <- ifelse(
    recall, dunif(params$recall, min = 0, max = 5), 1
  )

  log(prob_a) + log(prob_b) + log(prob_c)
    log(prob_pinvalid) + log(prob_recall)
}

log_prior_gamma <- function(params, mixture, recall) {

  prob_a <- dunif(params$alpha1, min = 0, max = 100)
  prob_b <-  dunif(params$beta1, min = 0, max = 100)
  out <- log(prob_a) + log(prob_b)
  out <- ifelse(mixture, out + dbeta(params$pinvalid, shape1 = 4, shape2 = 10), out)
  out <- ifelse(recall, out + dunif(params$recall, min = 0, max = 5), out)
  out
}

log_likel <- function(samples, mixture, recall,
                      model = c("nf", "beta", "skew_normal", "gamma"),
                      a_priors = list(mean = 4, sd = 1),
                      b_priors = list(mean = 1, sd = 0.5)) {
  map_dbl(
    seq_along(samples[[1]]), function(index) {
      params <- map(samples, ~ .[[index]])
      if (model == "nf") {
        log_prior <- log_prior_nf(params, mixture, recall, a_priors, b_priors)
      } else if (model == "skew_normal") {
        log_prior <- log_prior_skew_normal(params, mixture, recall)
      } else if (model == "gamma") {
        log_prior <- log_prior_gamma(params, mixture, recall)
      } else {
        log_prior <- log_prior_beta(params, mixture, recall)
      }
      params$lp_ - log_prior
    }
  )
}

## calculating DIC according to Gelman 2004 + Wiki
DIC_alt <- function(loglik) {
  D <- -2 * loglik
  Dbar <- mean(D)
  pD <- 0.5 * var(D)
  Dbar + pD
}
