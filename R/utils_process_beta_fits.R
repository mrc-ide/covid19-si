sample_inf_profile <- function(n, alpha1, beta1, max_shed, offset) {
  inf_times <- rbeta(n = n, shape1 = alpha1, shape2 = beta1)
  ((max_shed - offset) * inf_times) + offset
}

sample_si <- function(n = 1e4,
                      inf_samples,
                      shape_inc = params_real$inc_par2[["shape"]],
                      scale_inc = params_real$inc_par2[["scale"]]
                      ) {
  if (! is.matrix(inf_samples)) {
    inf_samples <- matrix(
      inf_samples, nrow = n, ncol = 1
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

  inc_mat + inf_samples
}

process_beta_fit <- function(pars, n = 1e4, samples,
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
    nrow = n, ncol = length(samples$alpha1)
  )
  for (i in seq_len(length(samples$alpha1))) {
    TOST_post[,i] <- sample_inf_profile(
    n, samples$alpha1[i], samples$beta1[i], max_shed, offset
   )
  }

  ## Serial Interval
  SI_best <- sample_si(n, TOST_bestpars)
  SI_mean <- sample_si(n, TOST_meanpars)
  SI_post <- sample_si(n, TOST_post)

  list(
    TOST_bestpars = TOST_bestpars,
    TOST_meanpars = TOST_meanpars,
    TOST_post = TOST_post,
    SI_bestpars = data.frame(
      SI = SI_best[, 1], TOST = TOST_bestpars
    ),
    SI_meanpars = data.frame(
      SI = SI_mean[, 1], TOST = TOST_meanpars
    ),
    SI_post = SI_post
  )
}


cowling_data <- readRDS("data/cowling_data_clean.rds") %>%
  mutate(si = as.numeric(si))%>%
   dplyr::rename(nu = onset_first_iso)%>%
   dplyr::filter(!is.na(nu))



### S4 Mixture
fit4mix <- readRDS("stanfits/release/scenario4a_mixture_beta.rds")
tab1_s4mix <- table1_fun(fit4mix)
samples_s4 <- process_beta_fit(
  tab1_s4mix, 1e4, rstan::extract(fit4mix), 21, -11
)
tab2_s4mix <- table2_fun(samples_s4)


TOST4 <- samples_s4$TOST_bestpars
p <- TOST_fig_fun(TOST4)
SI4 <- samples_s4$SI_bestpars$SI


SI4_con <- expected_SI_fun(
  SI = samples_s4$SI_bestpars,
  data = cowling_data,
  mixture = TRUE,
  recall = FALSE,
  isol = TRUE,
  tab1 = tab1_s4mix,
  n = 1e4
)

p <- SIcomp_fig_fun(SI1 = SI4, SI2 = SI4_con, data = cowling_data)
ggsave("figures/s4mixture_cowling_beta.png", p)


### S4 Mixture and recall
fit4recall <- readRDS(
  "stanfits/release/scenario4arecall_mixture_beta.rds"
)
tab1_s4recall <- table1_fun(fit4recall)

samples_s4recall <- process_beta_fit(
  tab1_s4recall, 1e4, rstan::extract(fit4recall), 21, -11
)
tab2_s4recall <- table2_fun(samples_s4recall)


TOST4 <- samples_s4recall$TOST_bestpars
p <- TOST_fig_fun(TOST4)
SI4recall <- samples_s4recall$SI_bestpars$SI


SI4recall_con <- expected_SI_fun(
  SI = samples_s4recall$SI_bestpars,
  data = cowling_data,
  mixture = TRUE,
  recall = TRUE,
  isol = TRUE,
  tab1 = tab1_s4recall,
  n = 1e4
)

p <- SIcomp_fig_fun(SI1 = SI4recall, SI2 = SI4recall_con, data = cowling_data)
ggsave("figures/s4recall_cowling_beta.png", p)
