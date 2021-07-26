## Returns TOST for
## (a) parameters with maximum posterior likelihood
## (b) at the mean of the posterior distribution
## (c) a distribution of distributions of TOST ; 1
## distribution for each parameter in the posterior
## distributions of xi, omega, and alpha  sampled jointly.
## tab1 is the output from fitted_params function
estimated_TOST_skew_normal <- function(tab1, n = 1e5, fit) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- sn::rsn(
    n = n, xi = best["a"], omega = best["b"], alpha = best["c"]
  )

  mu_params <- tab1$mean
  names(mu_params) <- rownames(tab1)
  TOST_meanpars <- sn::rsn(
    n = n, xi = mu_params["a"], omega = mu_params["b"],
    alpha = mu_params["c"]
  )

  TOST_post <- matrix(nrow = n, ncol = length(fit$a))
  for (i in 1:(length(fit$a))) {
    TOST_post[, i] <- sn::rsn(
      n = n, xi = fit$a[i], omega = fit$b[i], alpha = fit$c[i]
    )
  }

  TOST_post <- as.data.frame(TOST_post)

  list(
    TOST_bestpars = TOST_bestpars,
    TOST_meanpars = TOST_meanpars,
    TOST_post = TOST_post
  )
}
