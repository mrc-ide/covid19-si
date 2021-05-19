sample_inf_profile <- function(n, alpha1, beta1, max_shed, offset) {
  inf_times <- rbeta(n = n, shape1 = alpha1, shape2 = beta1)
  ((max_shed - offset) * inf_times) + offset
}

estimated_TOST_beta <- function(tab1, max_shed, offset, nsim, samples) {
    best <- tab1$best
    names(best) <- rownames(tab1)
    TOST_bestpars <- sample_inf_profile(
      nsim, best[["alpha1"]], best[["beta1"]], max_shed, offset
    )
    mean_pars <- tab1$mean
    names(mean_pars) <- rownames(tab1)
    TOST_meanpars <- sample_inf_profile(
      nsim, mean_pars[["alpha1"]], mean_pars[["beta1"]], max_shed, offset
    )

    ## Posterior distribution
    TOST_post <- matrix(
      nrow = nsim, ncol = length(samples$alpha1)
    )
    for (i in seq_len(length(samples$alpha1))) {
      TOST_post[,i] <- sample_inf_profile(
        nsim, samples$alpha1[i], samples$beta1[i], max_shed, offset
      )
    }

    TOST_post <- as.data.frame(TOST_post)

    list(
      TOST_bestpars = TOST_bestpars,
      TOST_meanpars = TOST_meanpars,
      TOST_post = TOST_post
    )
}
