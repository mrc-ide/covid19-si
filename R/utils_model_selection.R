## params is a named list of estimated parameters
log_prior_nf <- function(params, mixture, recall) {

  prob_a <- dnorm(params$a, mean = 4.28, sd = 0.74) / pnorm(0, mean = 4.28, sd = 0.74)
  prob_b <-  dnorm(params$b, mean = 1.44, sd = 0.12) / pnorm(0, mean = 1.44, sd = 0.12)
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


log_likel <- function(samples, mixture, recall, model = c("nf", "beta")) {
  map_dbl(
    seq_along(samples[[1]]), function(index) {
      params <- map(samples, ~ .[[index]])
      if (model == "nf") {
        log_prior <- log_prior_nf(params, mixture, recall)
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
