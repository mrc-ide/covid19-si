# extract fitted parameters
fitted_params <- function(fit, digits = 2) {
  tab1 <- as.data.frame(rstan::summary(fit)$summary)
  fit <- rstan::extract(fit)

  best_idx <- which.max(fit[["lp__"]])
  best <- unlist(map(fit, function(x) x[[best_idx]]))
  tab1 <- add_column(tab1, best = best, .after = 1)

  round(tab1, digits)
}


unconditional_si <- function(inf_samples,
                             shape_inc = params_real$inc_par2[["shape"]],
                             scale_inc = params_real$inc_par2[["scale"]]
                             ) {
  inc_samples <- rgamma(
    shape = shape_inc,
    scale = scale_inc,
    n = length(inf_samples)
  )
  list(
    inf = inf_samples,
    si = inc_samples + inf_samples
  )
}

## Return a named list where the names are the nus
## in the observe data set. This is so that
## (a) we can use s3 + recall etc i.e. make use of
## nu even when we don't condition on isolation,
## and (b) we have a consistent output format
## across all functions
unconditional_si_all <- function(obs, inf_samples) {

  x <- unconditional_si(inf_samples)
  nu_freq <- janitor::tabyl(obs$nu)
  colnames(nu_freq)[1] <- "nu"
  si_given_nu <- map(nu_freq$nu, function(nu) x)
  names(si_given_nu) <- nu_freq$nu
  si_given_nu
}

## un_si is the output of unconditional_si
##
conditional_si <- function(un_si, nu, nsim) {
  inf_samples <- un_si[[1]]
  idx <- which(inf_samples <= nu)
  filtered <- inf_samples[idx]
  filtered_si <- un_si[[2]][idx]

  if (length(filtered) == 0) return(NULL)
  inf_samples <- sample(filtered, nsim, replace = TRUE)
  si <- sample(filtered_si, nsim, replace = TRUE)
  list(inf = inf_samples, si = si)
}

## obs is observed data with column nu
## processed_fit is the output of process_beta_fit
conditional_si_all <- function(un_si, nsim) {
  si_given_nu <- imap(
    un_si, function(si, nu) {
      conditional_si(si, as.numeric(nu), nsim = nsim)
    }
  )
  keep(si_given_nu, function(x) !is.null(x[["si"]]))
}

## si_given_nu is the output of beta_fit_si_all, a list of conditional
## SIs for each nu in the data
## Returns a named list
conditional_si_pooled <- function(obs, si_given_nu, nsim) {
  ## Under some -ve nus, we have 0 possible SIs, we wil filter them out
  ## here.
  si_given_nu <- map(si_given_nu, function(x) x[["si"]])

  ## Renormalise weights over this smaller set of nus
  nu_freq <- janitor::tabyl(obs$nu[obs$nu %in% as.numeric(names(si_given_nu))])
  ## Number of samples to be drawn is nsim * percent
  nu_freq$nu_weights <- ceiling(
    nu_freq$percent * nsim
  )
  ## This can give you more or less than nsim
  ## SIs in all. You can resample when pooling.
  map2(
    si_given_nu, nu_freq$nu_weights, function(x, size) {
      sample(x, size, replace = TRUE)
    }
  )
}

leaky_si <- function(un_si, nu, pleak, nsim) {
  inf_samples <- un_si[[1]]

  ## Split inf_times into before and
  ## after nu. And then from the before
  ## group, select 1 - pleak samples
  ## and from the after group, select
  ## pleak samples
  before_nu <- which(inf_samples <= nu)
  after_nu <- which(inf_samples > nu)
  ## Do it this way rather than sampling the vector
  ## directly so that we can use the index to select
  ## SIs
  leaky <- sample(
    length(after_nu), size = floor(nsim * pleak),
    replace = TRUE
  )
  not_leaky <- sample(
    length(before_nu), size = floor(nsim * (1 - pleak)),
    replace = TRUE
  )
  filtered <- c(
    inf_samples[before_nu][not_leaky],
    inf_samples[after_nu][leaky]
  )
  filtered_si <- c(
    un_si[[2]][before_nu][not_leaky],
    un_si[[2]][after_nu][leaky]
  )
  list(
    inf = filtered, si = filtered_si,
    leaky_inf = after_nu[leaky],
    not_leaky_inf = before_nu[not_leaky],
    leaky_si = un_si[[2]][after_nu][leaky],
    not_leaky_si = un_si[[2]][before_nu][not_leaky]
  )
}

leaky_si_all <- function(un_si, pleak, nsim) {
  si_given_nu <- imap(
    un_si, function(si, nu) {
      leaky_si(si, as.numeric(nu), pleak, nsim = nsim)
    }
  )
  keep(si_given_nu, function(x) !is.null(x[["si"]]))
}


## predicted observed SIs (under assumed biases)
## currently assumes recall and isolation biases only affect valid SIs
##
## This function will then apply relevant biases to
estimated_SI <- function(obs, inf_times, mixture, recall,
                         isol, leaky = FALSE, tab1, nsim = 1e4, tmin = -20) {

  if (leaky & isol) {
    stop("leaky and isolation should both not be true")
  }


  ## First simulate unconditional SI and then apply biases
  un_si <- unconditional_si_all(obs, inf_times)
  pleak_par <- ifelse(leaky, tab1["pleak", "best"], 0)
  leaky <- leaky_si_all(un_si, pleak_par, nsim)
  # with isolation bias
  if (isol) {
    with_iso <- conditional_si_all(leaky, nsim)
  } else {
    with_iso <- leaky
  }
  with_iso <- map(with_iso, ~ .[["si"]])
  ## with recall bias
  recall_par <- ifelse(recall, tab1["recall", "best"], 0)

  with_recall <- imap(
    with_iso, function(si, nu) {
      nu <- as.numeric(nu)
      precall <- exp(-recall_par * abs(si - nu))
      ## Keep the same list structure
      list(
        si = sample(si, size = nsim, replace = TRUE, prob = precall)
      )
    }
  )
  ## Now pool, this is still a named list with names being nu values
  un_si <- conditional_si_pooled(obs, un_si, nsim)
  pooled_unc <- unname(unlist(un_si))

  si <- conditional_si_pooled(obs, with_recall, nsim)
  pooled_si <- unname(unlist(si))

  ## Mixture
  pinvalid <- ifelse(mixture, tab1["pinvalid", "best"], 0)
  toss <- runif(nsim)
  valid <- which(toss >= pinvalid)
  invalid_si <- runif(nsim, tmin, 40)

  mixed <- c(pooled_si[valid], invalid_si[-valid])

  list(
    unconditional = pooled_unc, conditional = mixed
  )
}
