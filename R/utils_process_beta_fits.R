sample_inf_profile <- function(n, alpha1, beta1, max_shed, offset) {
  inf_times <- rbeta(n = n, shape1 = alpha1, shape2 = beta1)
  ((max_shed - offset) * inf_times) + offset
}

process_beta_fit <- function(pars, n = 1e5, samples,
                             shape_inc = params_real$inc_par2[["shape"]],
                             scale_inc = params_real$inc_par2[["scale"]],
                             max_shed, offset) {

  best <- pars$best
  names(best) <- rownames(pars)

  TOST_bestpars <- sample_inf_profile(
    n, best[["alpha1"]], best[["beta1"]], max_shed, offset
  )

  mean_pars <- pars$mean
  names(mean_pars) <- rownames(pars)
  TOST_meanpars <- sample_inf_profile(
    n, mean_pars[["alpha1"]], mean_pars[["beta1"]], max_shed, offset
  )

  ## Posterior distribution
  TOST_post <- matrix(
    ncol = length(samples$alpha1), nrow = n
  )
  for (i in seq_len(length(samples$alpha1))) {
    TOST_post[,i] <- sample_inf_profile(
    n, samples$alpha1[i], samples$beta1[i], max_shed, offset
    )
  }

  TOST_post <- as.data.frame(TOST_post)

  ## Serial Interval
  inc_mat <- as.data.frame(matrix(rgamma(shape = shape_inc,
                                         scale = scale_inc,
                                         n = n*length(fit$a)),
                                  nrow = dim(TOST_post)[1],
                                  ncol = dim(TOST_post)[2]))
  SI_best <- inc_mat[,1] + TOST_bestpars
  SI_mean <- inc_mat[,1] + TOST_meanpars
  SI_post <- TOST_post + inc_mat

  list(
    TOST_bestpars = TOST_bestpars,
    TOST_meanpars = TOST_meanpars,
    TOST_post = TOST_post,
    SI_bestpars = data.frame(SI = SI_best, TOST = TOST_bestpars),
    SI_meanpars = data.frame(SI = SI_mean, TOST = TOST_meanpars),
    SI_post = SI_post
  )
}


fit4mix <- readRDS("stanfits/release/scenario4a_mixture_beta.rds")
fit4recall <- readRDS("stanfits/release/scenario4arecall_mixture_beta.rds")
tab1_s4mix <- table1_fun(fit4mix)
tab1_s4recall <- table1_fun(fit4recall)
