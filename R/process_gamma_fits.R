## Returns TOST for
## (a) parameters with maximum posterior likelihood
## (b) at the mean of the posterior distribution
## (c) a distribution of distributions of TOST ; 1
## distribution for each parameter in the posterior
## distributions of xi, omega, and alpha  sampled jointly.
## tab1 is the output from fitted_params function
estimated_TOST_gamma <- function(tab1, n = 1e4, fit, offset = 20) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- rgamma(
    n = n, shape = best[["alpha1"]], rate = best[["beta1"]]
  ) - offset

  mu_params <- tab1$mean
  names(mu_params) <- rownames(tab1)
  TOST_meanpars <- rgamma(
    n = n, shape = mu_params["alpha1"], rate = mu_params["beta1"],
  ) - offset

  TOST_post <- matrix(nrow = n, ncol = length(fit$alpha1))
  for (i in 1:(length(fit$alpha1))) {
    TOST_post[, i] <- rgamma(
      n = n, shape = fit$alpha1[i], rate = fit$beta1[i]
    ) - offset
  }

  TOST_post <- as.data.frame(TOST_post)

  list(
    TOST_bestpars = TOST_bestpars,
    TOST_meanpars = TOST_meanpars,
    TOST_post = TOST_post
  )
}

check <- "\U2713"
meta_model <- "s3s4"
fit_dir <- "stanfits/gamma/s3s4"
outdir <- "processed_stanfits/gamma/s3s4"
figs_dir <- "figures/gamma/s3s4"

if (grepl("discrete_pairs", meta_model)) {
  obs_data <- data_discrete_pairs
} else if (grepl("s3s4", meta_model)) {
  obs_data <- data_s3_s4mix
} else  {
  obs_data <- cowling_data
}

##obs_data <- data_discrete_pairs_s3_s4mix

if (! dir.exists("processed_stanfits")) dir.create("processed_stanfits")
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (! dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

infiles <- glue("{fit_dir}/{model_features$model_prefix}_gamma_fit.rds")
index <- map_lgl(infiles, file.exists)
model_features <- model_features[index, ]
infiles <- infiles[index]

fits <- map(infiles, readRDS)

names(fits) <- model_features$model_prefix

table1 <- map2(
  fits, model_features$model_prefix,
  function(fit, model_prefix) {
    tab1 <- fitted_params(fit)
    message(
      glue("{check} Extracted parameters for model {model_prefix}")
    )
    saveRDS(
      tab1, glue("{outdir}/{model_prefix}_gamma_tab1.rds")
    )
    tab1
  }
)

## table1 <- map(
##   model_features$model_prefix,
##   function(x) readRDS(glue("{outdir}/{x}_gamma_tab1.rds"))
## )

samples_tost <- map2(
  table1, fits, function(tab1, fit) {
    estimated_TOST_gamma(
      tab1, n = 1e4, rstan::extract(fit)
    )
  }
)

saveRDS(samples_tost, glue("{outdir}/samples_tost_gamma.rds"))

best_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1, model_prefix = model_features$model_prefix
  ), function(tost, mixture, recall, isol, tab1, model_prefix) {
    message(
      glue("{check} Sampling best SI for model {model_prefix}")
    )
    estimated_SI(
      obs_data, tost$TOST_bestpars, mixture = mixture,
      recall = recall, isol = isol, leaky = FALSE, tab1 = tab1
    )
  }
)

saveRDS(best_si, glue("{outdir}/best_si_gamma.rds"))

mean_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1, model_prefix = model_features$model_prefix
    ), function(tost, mixture, recall, isol, tab1, model_prefix) {
    message(
      glue("{check} Sampling mean SI for model {model_prefix}")
    )
    estimated_SI(
      obs_data, tost$TOST_meanpars, mixture = mixture,
      recall = recall, isol = isol, tab1 = tab1
    )
  }
)

saveRDS(mean_si, glue("{outdir}/mean_si_gamma.rds"))

post_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1, model_prefix = model_features$model_prefix
  ), function(tost, mixture, recall, isol, tab1, model_prefix) {
    message(
      glue("{check} Posterior SI distribution for model {model_prefix}")
    )
    ## Posterior has many samples here. Thin it down.
    ## Each column is a sample of length 10000 and there are 10000
    ## of them
    cols <- seq(1, ncol(tost$TOST_post), 2)
    tost$TOST_post <- tost$TOST_post[, cols]
    x <- apply(
      tost$TOST_post, 2, function(inf_samples) {
        estimated_SI(
          obs_data, inf_samples, mixture = mixture,
          recall = recall, isol = isol, tab1 = tab1
        )
      }
    )
    ## Extract conditional SI and make a matrix/data.frame
    map(x, ~ .[["unconditional"]]) %>% do.call(what = 'rbind')
  }
)

saveRDS(post_si, glue("{outdir}/gamma_post_si.rds"))
## post_si <- readRDS("{outdir}/gamma_post_si.rds")

table2 <- pmap(
  list(
    mean_si = mean_si, best_si = best_si, post = post_si,
    tost = samples_tost,
    model_prefix = model_features$model_prefix
  ),
  function(mean_si, best_si, post, tost, model_prefix) {
    message(
      glue("{check} Constructing table 2 for model {model_prefix}")
    )
    samples_si <- list(
      SI_meanpars = list(SI = mean_si[["unconditional"]]),
      SI_bestpars = list(SI = best_si[["unconditional"]]),
      SI_post = post
    )
    ## table 2 - summary stats for sampled distributions
    tab2 <- tost_si_summary(tost, samples_si)
    saveRDS(
      tab2, glue("{outdir}/{model_prefix}_gamma_tab2.rds")
    )
    tab2
  }
)

walk2(
  samples_tost, model_features$model_prefix,
  function(tost, model_prefix) {
  p1 <- TOST_figure(tost$TOST_bestpars)
  ggsave(
    filename = glue("{figs_dir}/{model_prefix}_gamma_tost.png"), p1
  )
})


walk2(
  best_si, model_features$model_prefix,
  function(si, model_prefix) {
  psi <- SI_figure(si[[2]], obs_data)
  ggsave(
    filename = glue("{figs_dir}/{model_prefix}_gamma_si.png"), psi
  )
})

x <- as.list(model_features)
x <- append(x, values = list(si = best_si))
pwalk(
  x, function(mixture, recall,
              right_bias, model_prefix, si) {
    ## If none of these are true, then conditional
    ## and unconditional SIs would be the same.
    message("model_prefix ", model_prefix)
    if (mixture |  recall | right_bias) {
        psi2 <- plot_both_SIs(
          SI1 = si[[1]], SI2 = si[[2]],
          data = obs_data
        )
        ggsave(
          filename = glue("{figs_dir}/{model_prefix}_gamma_si.png"), psi2
        )
    }
  }
)

x <- as.list(model_features)
x <- append(x, values = list(fits = fits))

dic <- pmap_dbl(
  x, function(mixture, recall, right_bias, model_prefix, fits) {
    samples <- rstan::extract(fits)
    DIC_alt(
      log_likel(samples, mixture, recall, "gamma")
    )
  }
)
names(dic) <- x$model_prefix
## To be executed after getting Table 2s from
## various models
overall_table2 <- map_dfr(
  model_features$model_prefix,
  function(model_prefix) {
    out <- readRDS(glue("{outdir}/{model_prefix}_gamma_tab2.rds"))
    out <- tibble::rownames_to_column(out, var = "param")
    out$DIC <- dic[[model_prefix]]
    out$model <- model_prefix

    out
  }
)

overall_table2$formatted_pars <-
  glue(
    "{overall_table2$best_pars} ",
    "({overall_table2$CrI_2.5} - {overall_table2$CrI_97.5})"
  )

x <- overall_table2[ ,c("param", "formatted_pars", "model", "DIC")]
x <- tidyr::spread(x, key = param, value = formatted_pars)
x <- arrange(x, dic)
saveRDS(x, glue("{outdir}/gamma_overall_table2.rds"))



## For manuscript
## for_ms <- select(x, model, mean_igamma, sd_inf, mean_si, sd_si, DIC)
## for_ms <- arrange(for_ms, DIC)
## for_ms <- left_join(
##   for_ms, model_features, by = c("model" = "model_prefix")
## )
## for_ms <- select(
##  for_ms, mixture, recall, right_bias, mean_inf, sd_inf, mean_si, sd_si, DIC
## )
## stargazer::stargazer(for_ms, summary = FALSE, rownames = FALSE)
