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
                             scale_inc = params_real$inc_par2[["scale"]],
                             nsim = 1e4
                             ) {
  if (! is.matrix(inf_samples)) {
    inf_samples <- matrix(
      inf_samples, nrow = nsim, ncol = 1
    )
  }
  inc_mat <- matrix(
      rgamma(
        shape = shape_inc,
        scale = scale_inc,
        n = nrow(inf_samples) *
          ncol(inf_samples)
      ),
      nrow = nrow(inf_samples),
      ncol = ncol(inf_samples)
    )

  list(inf = inf_samples, si = inc_mat + inf_samples)
}

## Simulate serial interval distribution given nu
conditional_si <- function(nu, inf_samples,
                           shape_inc = params_real$inc_par2[["shape"]],
                           scale_inc = params_real$inc_par2[["scale"]],
                           nsim) {

  filtered <- inf_samples[inf_samples <= nu]
  if (length(filtered) == 0) return(NULL)
  inf_samples <- sample(filtered, nsim, replace = TRUE)
  inc_times <-  rgamma(n = nsim, shape = shape_inc, scale = scale_inc)
  list(inf = inf_samples, si = inc_mat + inf_samples)
}

## obs is observed data with column nu
## processed_fit is the output of process_beta_fit
conditional_si_all <- function(obs, inf_samples, nsim) {
  nu_freq <- janitor::tabyl(obs$nu)
  colnames(nu_freq)[1] <- "nu"
  si_given_nu <- map(
    nu_freq$nu, function(nu) {
      conditional_si(
        nu = nu, inf_samples = inf_samples, nsim = nsim
      )
    }
  )
  names(si_given_nu) <- nu_freq$nu
  si_given_nu
}
## si_given_nu is the output of beta_fit_si_all, a list of conditional
## SIs for each nu in the data
conditional_si_pooled <- function(obs, si_given_nu) {
  ## Under some -ve nus, we have 0 possible SIs, we wil filter them out
  ## here.
  si_given_nu <- keep(si_given_nu, function(x) !is.null(x))
  ## Renormalise weights over this smaller set of nus
  nu_freq <- janitor::tabyl(obs$nu[obs$nu %in% as.numeric(names(si_given_nu))])
  ## Now sample from this list.
  nu_weights <- sample(
    names(si_given_nu), nsim, replace = TRUE, prob = nu_freq$percent
  )
  nu_weights <- janitor::tabyl(nu_weights)
  si <- map2(
    si_given_nu, nu_weights$n, function(x, size) {
      sample(x, size)
    }
  )
  unname(unlist(si))
}

