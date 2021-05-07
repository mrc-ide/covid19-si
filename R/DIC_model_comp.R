# model comparison using DIC
rstan::expose_stan_functions("stan-models/likelihoods_nf.stan")

# import fits
fit3 <- rstan::extract(readRDS("stanfits/release/scenario3a_mixture_nf.rds"))
fit3_recall <- rstan::extract(readRDS("stanfits/release/scenario3arecall_mixture_nf.rds"))
fit4 <- rstan::extract(readRDS("stanfits/release/scenario4a_mixture_nf.rds"))
fit4_recall <- rstan::extract(readRDS("stanfits/release/scenario4arecall_mixture_nf.rds"))
fit3_nomix_recall <- rstan::extract(readRDS("stanfits/release/scenario3arecall_no_mixture_nf.rds"))
fit4_nomix_recall <- rstan::extract(readRDS("stanfits/release/scenario4arecall_no_mixture_nf.rds"))
fit3_leftbias <- rstan::extract(readRDS("stanfits/release/scenario3a_mixture_left_bias_nf.rds"))
fit3_nomix_leftbias <- rstan::extract(readRDS("stanfits/release/scenario3a_nomixture_left_bias_nf.stan.rds"))

## calculating log(prior)

log_prior_func_basic <- function(a, b, c, tmax){
  prob_a <- dnorm(a, mean = 4.28, sd = 0.74) / pnorm(0, mean = 4.28, sd = 0.74)
  prob_b <-  dnorm(b, mean = 1.44, sd = 0.12) / pnorm(0, mean = 1.44, sd = 0.12)
  prob_c <- dunif(c, min = 0, max = 1)
  prob_tmax <- dunif(tmax, min = -20, max = 10)
  
  log_prior <- log(prob_a) + log(prob_b) + log(prob_c) + log(prob_tmax)
  
}

log_prior_func_mix <- function(a, b, c, tmax, pinvalid){
  prob_a <- dnorm(a, mean = 4.28, sd = 0.74) / pnorm(0, mean = 4.28, sd = 0.74)
  prob_b <-  dnorm(b, mean = 1.44, sd = 0.12) / pnorm(0, mean = 1.44, sd = 0.12)
  prob_c <- dunif(c, min = 0, max = 1)
  prob_tmax <- dunif(tmax, min = -20, max = 10)
  prob_pinvalid <- dbeta(pinvalid, shape1 = 4, shape2 = 10)
  
  log_prior <- log(prob_a) + log(prob_b) + log(prob_c) + log(prob_tmax) + log(prob_pinvalid)
}

log_prior_func_recall <- function(a, b, c, tmax, recall){
  prob_a <- dnorm(a, mean = 4.28, sd = 0.74) / pnorm(0, mean = 4.28, sd = 0.74)
  prob_b <-  dnorm(b, mean = 1.44, sd = 0.12) / pnorm(0, mean = 1.44, sd = 0.12)
  prob_c <- dunif(c, min = 0, max = 1)
  prob_tmax <- dunif(tmax, min = -20, max = 10)
  prob_recall <- dunif(recall, min = 0, max = 5)
  
  log_prior <- log(prob_a) + log(prob_b) + log(prob_c) + log(prob_tmax) + log(prob_recall)
  
}

log_prior_func_mix_recall <- function(a, b, c, tmax, pinvalid, recall){
  prob_a <- dnorm(a, mean = 4.28, sd = 0.74) / pnorm(0, mean = 4.28, sd = 0.74)
  prob_b <-  dnorm(b, mean = 1.44, sd = 0.12) / pnorm(0, mean = 1.44, sd = 0.12)
  prob_c <- dunif(c, min = 0, max = 1)
  prob_tmax <- dunif(tmax, min = -20, max = 10)
  prob_pinvalid <- dbeta(pinvalid, shape1 = 4, shape2 = 10)
  prob_recall <- dunif(recall, min = 0, max = 5)
  
  log_prior <- log(prob_a) + log(prob_b) + log(prob_c) + log(prob_tmax) + log(prob_pinvalid) + log(prob_recall)
  
}

# generating the log likelihood by subtracting log(prior) from lp__
loglik_func_basic <- function(fit){
  logprior <- vector(length = length(fit$a))
  for(i in 1:(length(fit$a))){
    logprior[i] <- log_prior_func_basic(a = fit$a[i],
                                      b = fit$b[i],
                                      c = fit$c[i],
                                      tmax = fit$tmax[i])
  }
  loglik <- fit$lp__ - logprior
}

loglik_func_mix <- function(fit){
  logprior <- vector(length = length(fit$a))
for(i in 1:(length(fit$a))){
 logprior[i] <- log_prior_func_mix(a = fit$a[i],
                               b = fit$b[i],
                               c = fit$c[i],
                               tmax = fit$tmax[i],
                               pinvalid = fit$pinvalid[i])
}
loglik <- fit$lp__ - logprior
}

loglik_func_recall <- function(fit){
  logprior <- vector(length = length(fit$a))
  for(i in 1:(length(fit$a))){
    logprior[i] <- log_prior_func_recall(a = fit$a[i],
                                             b = fit$b[i],
                                             c = fit$c[i],
                                             tmax = fit$tmax[i],
                                             recall = fit$recall[i])
  }
  loglik <- fit$lp__ - logprior
}


loglik_func_mix_recall <- function(fit){
  logprior <- vector(length = length(fit$a))
  for(i in 1:(length(fit$a))){
    logprior[i] <- log_prior_func_mix_recall(a = fit$a[i],
                                      b = fit$b[i],
                                      c = fit$c[i],
                                      tmax = fit$tmax[i],
                                      pinvalid = fit$pinvalid[i],
                                      recall = fit$recall[i])
  }
  loglik <- fit$lp__ - logprior
}

## calculating DIC according to Spiegelhalter 2002
DIC_func <- function (loglik, loglik_at_mean) {
  Dbar <- -2 * mean(loglik)
  pD <- Dbar - (-2 * loglik_at_mean)
  DIC <- Dbar + pD
  return(DIC)
}

## calculating DIC according to Gelman 2004 + Wiki
DIC_alt_func <- function(loglik) {
  D <- -2 * loglik
  Dbar <- mean(D)
  pD <- 0.5 * var(D)
  DIC <- Dbar + pD
  return(DIC)
}


loglik3 <- loglik_func_mix(fit3)
loglik3_recall <- loglik_func_mix_recall(fit3_recall)
loglik3_nomix_recall <- loglik_func_recall(fit3_nomix_recall)
loglik4 <- loglik_func_mix(fit4)
loglik4_recall <- loglik_func_mix_recall(fit4_recall)
loglik4_nomix_recall <- loglik_func_recall(fit4_nomix_recall)
loglik3_leftbias <- loglik_func_mix(fit3_leftbias)
loglik3_nomix_leftbias <- loglik_func_basic(fit3_nomix_leftbias)

DIC3 <- DIC_alt_func(loglik3)
DIC3_recall <- DIC_alt_func(loglik3_recall)
DIC3_nomix_recall <- DIC_alt_func(loglik3_nomix_recall)
DIC4 <- DIC_alt_func(loglik4)
DIC4_recall <- DIC_alt_func(loglik4_recall)
DIC4_nomix_recall <- DIC_alt_func(loglik4_nomix_recall)
DIC3_leftbias <- DIC_alt_func(loglik3_leftbias)
DIC3_nomix_leftbias <- DIC_alt_func(loglik3_nomix_leftbias)

