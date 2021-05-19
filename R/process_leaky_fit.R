fit <- readRDS("stanfits/scenario3a_leaky_nf_fit.rds")
tab1 <- fitted_params(fit)
tost <- estimated_TOST_nf(
  tab1, taus = seq(-20, 40, 0.1), n = 1e4,
  rstan::extract(fit)
)
best_si <-  estimated_SI(
  cowling_data, tost$TOST_bestpars, mixture = FALSE,
  recall = FALSE, isol = FALSE, leaky = TRUE, tab1 = tab1
)
mean_si <- estimated_SI(
  cowling_data, tost$TOST_meanpars, mixture = mixture,
  recall = recall, isol = isol, tab1 = tab1
)
x <- apply(
  tost$TOST_post, 2, function(inf_samples) {
    estimated_SI(
      cowling_data, inf_samples, mixture = FALSE,
      recall = FALSE, isol = FALSE, leaky = TRUE, tab1 = tab1
    )
  }
)
## Extract conditional SI and make a matrix/data.frame

post_si <- map(x, ~ .[[2]]) %>% do.call(what = 'rbind')
samples_si <- list(
  SI_meanpars = list(SI = mean_si[[2]]),
  SI_bestpars = list(SI = best_si[[2]]),
  SI_post = post
)
## table 2 - summary stats for sampled distributions
tab2 <- tost_si_summary(tost, samples_si)

psi <- SI_figure(best_si[[2]], cowling_data)

save_plot(
  filename = glue("figures/{model_prefix}_nf_si.pdf"), psi
)
