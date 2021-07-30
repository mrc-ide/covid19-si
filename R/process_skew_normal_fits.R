source("R/utils_process_skew_normal_fits.R")
check <- "\U2713"
meta_model <- "discrete_pairs"
fit_dir <- "stanfits/skew_normal/discrete_pairs"
outdir <- "processed_stanfits/skew_normal/discrete_pairs"
figs_dir <- "figures/skew_normal/releasediscrete_pairs"

if (grepl("discrete_pairs", meta_model)) {
  obs_data <- data_discrete_pairs
} else if (grepl("s3s4mix", meta_model)) {
   obs_data <- data_s3_s4mix
} else  {
  obs_data <- cowling_data
}

##obs_data <- data_discrete_pairs_s3_s4mix

if (! dir.exists("processed_stanfits")) dir.create("processed_stanfits")
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (! dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

infiles <- glue("{fit_dir}/{model_features$model_prefix}_skew_normal_fit.rds")
index <- map_lgl(infiles, file.exists)
model_features <- model_features[index, ]
infiles <- infiles[index]

fits <- map(infiles, readRDS)

table1 <- map2(
  fits, model_features$model_prefix,
  function(fit, model_prefix) {
    tab1 <- fitted_params(fit)
    message(
      glue("{check} Extracted parameters for model {model_prefix}")
    )
    saveRDS(
      tab1, glue("{outdir}/{model_prefix}_skew_normal_tab1.rds")
    )
    tab1
  }
)

## table1 <- map(
##   model_features$model_prefix,
##   function(x) readRDS(glue("{outdir}/{x}_skew_normal_tab1.rds"))
## )

samples_tost <- map2(
  table1, fits, function(tab1, fit) {
    estimated_TOST_skew_normal(
      tab1, n = 1e4, rstan::extract(fit)
    )
  }
)

saveRDS(samples_tost, glue("{outdir}/samples_tost_skew_normal.rds"))

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

saveRDS(best_si, glue("{outdir}/best_si_skew_normal.rds"))

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

saveRDS(mean_si, glue("{outdir}/mean_si_skew_normal.rds"))

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

saveRDS(post_si, glue("{outdir}/skew_normal_post_si.rds"))
## post_si <- readRDS("{outdir}/skew_normal_post_si.rds")

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
      tab2, glue("{outdir}/{model_prefix}_skew_normal_tab2.rds")
    )
    tab2
  }
)

walk2(
  samples_tost, model_features$model_prefix,
  function(tost, model_prefix) {
  p1 <- TOST_figure(tost$TOST_bestpars)
  ggsave(
    filename = glue("{figs_dir}/{model_prefix}_skew_normal_tost.png"), p1
  )
})


walk2(
  best_si, model_features$model_prefix,
  function(si, model_prefix) {
  psi <- SI_figure(si[[2]], obs_data)
  ggsave(
    filename = glue("{figs_dir}/{model_prefix}_skew_normal_si.png"), psi
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
          filename = glue("{figs_dir}/{model_prefix}_skew_normal_si.png"), psi2
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
      log_likel(samples, mixture, recall, "skew_normal")
    )
  }
)
names(dic) <- x$model_prefix
## To be executed after getting Table 2s from
## various models
overall_table2 <- map_dfr(
  model_features$model_prefix,
  function(model_prefix) {
    out <- readRDS(glue("{outdir}/{model_prefix}_skew_normal_tab2.rds"))
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
saveRDS(x, glue("{outdir}/skew_normal_overall_table2.rds"))



## For manuscript
## for_ms <- select(x, model, mean_iskew_normal, sd_inf, mean_si, sd_si, DIC)
## for_ms <- arrange(for_ms, DIC)
## for_ms <- left_join(
##   for_ms, model_features, by = c("model" = "model_prefix")
## )
## for_ms <- select(
##  for_ms, mixture, recall, right_bias, mean_inf, sd_inf, mean_si, sd_si, DIC
## )
## stargazer::stargazer(for_ms, summary = FALSE, rownames = FALSE)
