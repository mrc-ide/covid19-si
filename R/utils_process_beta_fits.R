sample_inf_profile <- function(n, alpha1, beta1, max_shed, offset) {
  inf_times <- rbeta(n = n, shape1 = alpha1, shape2 = beta1)
  ((max_shed - offset) * inf_times) + offset
}

unconditional_si <- function(n = 1e4,
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

## Simulate serial interval distribution given nu
conditional_si <- function(nu, inf_times,
                           shape_inc = params_real$inc_par2[["shape"]],
                           scale_inc = params_real$inc_par2[["scale"]],
                           nsim) {

  filtered <- inf_times[inf_times <= nu]
  if (length(filtered) == 0) return(NULL)
  inf_times <- sample(filtered, nsim, replace = TRUE)
  inc_times <-  rgamma(n = nsim, shape = shape_inc, scale = scale_inc)
  inf_times + inc_times
}

## obs is observed data with column nu
## processed_fit is the output of process_beta_fit
beta_fit_si_all <- function(obs, processed_fit) {
  nu_freq <- janitor::tabyl(obs$nu)
  colnames(nu_freq)[1] <- "nu"
  si_given_nu <- map(
    nu_freq$nu, function(nu) {
      conditional_si(
        nu = nu, inf_times = processed_fit$TOST_bestpars, nsim = 1e4
      )
    }
  )
  names(si_given_nu) <- nu_freq$nu
  si_given_nu
}
## si_given_nu is the output of beta_fit_si_all, a list of conditional
## SIs for each nu in the data
beta_fit_si_pooled <- function(obs, si_given_nu) {
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
  SI_best <- unconditional_si(n, TOST_bestpars)
  SI_mean <- unconditional_si(n, TOST_meanpars)
  SI_post <- unconditional_si(n, TOST_post)

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


### S3 Mixture
fit3mix <- readRDS("stanfits/release/scenario3a_mixture_beta.rds")
tab1_s3mix <- table1_fun(fit3mix)
samples_s3 <- process_beta_fit(
  tab1_s3mix, 1e4, rstan::extract(fit3mix), 21, -3
)
tab2_s3mix <- table2_fun(samples_s3)
TOST3 <- samples_s3$TOST_bestpars
p <- TOST_fig_fun(TOST3)
SI3 <- samples_s3$SI_bestpars$SI
p <- SI_fig_fun(SI3, cowling_data)
cowplot::save_plot("figures/scenario3a_mixture_beta_si.png", p)




### S4 Mixture
fit4mix <- readRDS("stanfits/release/scenario4a_mixture_beta.rds")
tab1_s4mix <- table1_fun(fit4mix)
samples_s4 <- process_beta_fit(
  tab1_s4mix, 1e4, rstan::extract(fit4mix), 21, -11
)
tab2_s4mix <- table2_fun(samples_s4)
TOST4 <- samples_s4$TOST_bestpars
p <- TOST_fig_fun(TOST4)

s4_si_all <- beta_fit_si_all(cowling_data, samples_s4)
s3_si_qntls <- map_dfr(si_all, quantile, .id = "nu")
s3_si_qntls$nu <- as.numeric(s3_si_qntls$nu)


psi_given_nu <- ggplot(s3_si_qntls) +
  geom_linerange(aes(y = nu, xmin = `25%`, xmax = `75%`)) +
  geom_point(aes(y = nu, x = `50%`)) +
  xlab("Serial Interval (Median and IQR)") +
  ylab("Time from symptom onset to isolation") +
  theme_minimal()

ggsave("figures/s4mixture_si_given_nu_beta.png", psi_given_nu)

s4_si_pooled <- beta_fit_si_pooled(cowling_data, s4_si_all)

p <- SIcomp_fig_fun(
  SI1 = samples_s4$SI_bestpar$SI, SI2 = s4_si_pooled,
  data = cowling_data
)

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
