model_features <- list(
  "mixture" = c(TRUE, FALSE),
  ##"left_bias" = c(TRUE, FALSE),
  "recall"  = c(TRUE, FALSE),
  "right_bias" = c(TRUE, FALSE)
)
model_features <- expand.grid(model_features)
model_features$model_prefix <- ifelse(
  model_features$`right_bias`, "scenario4a", "scenario3a"
)

model_features$model_prefix <-ifelse(
  model_features$mixture,
  glue::glue("{model_features$model_prefix}_mixture"),
  model_features$model_prefix
)

## model_features$model_prefix <-ifelse(
##   model_features$`left_bias`,
##   glue::glue("{model_features$model_prefix}_leftbias"),
##   model_features$model_prefix
## )

model_features$model_prefix <-ifelse(
  model_features$`recall`,
  glue::glue("{model_features$model_prefix}_recall"),
  model_features$model_prefix
)

## Need to rename fit files, at the moment only 1 has the right
## name
##index <- c(12, 15, 16)
index <- map_lgl(
  model_features$model_prefix,
  function(model_prefix) file.exists(glue("stanfits/{model_prefix}_beta_fit.rds"))
)

model_features <- model_features[index, ]

if (! dir.exists("processed_stanfits")) dir.create("processed_stanfits")
check <- "\U2713"

fits <- map(
  model_features$model_prefix,
  function(model_prefix) {
    message("Reading fit file for ", model_prefix)
    readRDS(glue("stanfits/{model_prefix}_beta_fit.rds"))
  }
)

table1 <- map2(
  fits, model_features$model_prefix,
  function(fit, model_prefix) {
    tab1 <- fitted_params(fit)
    message(
      glue("{check} Extracted parameters for model {model_prefix}")
    )
    saveRDS(
      tab1, glue("processed_stanfits/{model_prefix}_beta_tab1.rds")
    )
    tab1
  }
)


samples_tost <- map2(
  table1, fits, function(tab1, fit) {
    estimated_TOST_beta(
      tab1, max_shed = 21, offset = -11, nsim = 1e4,
      rstan::extract(fit)
    )
  }
)

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
      cowling_data, tost$TOST_bestpars, mixture = mixture,
      recall = recall, isol = isol, tab1 = tab1
    )
  }
)

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
      cowling_data, tost$TOST_meanpars, mixture = mixture,
      recall = recall, isol = isol, tab1 = tab1
    )
  }
)

post_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1, model_prefix = model_features$model_prefix
  ), function(tost, mixture, recall, isol, tab1, model_prefix) {
    message(
      glue("{check} Posterior SI distribution for model {model_prefix}")
    )
    x <- apply(
      tost$TOST_post, 2, function(ibeta_samples) {
        estimated_SI(
          cowling_data, ibeta_samples, mixture = mixture,
          recall = recall, isol = isol, tab1 = tab1
        )
      }
    )
    ## Extract conditional SI and make a matrix/data.frame
    map(x, ~ .[[2]]) %>% do.call(what = 'rbind')
  }
)

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
        SI_meanpars = list(SI = mean_si[[2]]),
        SI_bestpars = list(SI = best_si[[2]]),
        SI_post = post
      )
      ## table 2 - summary stats for sampled distributions
      tab2 <- tost_si_summary(tost, samples_si)
      saveRDS(
        tab2, glue("processed_stanfits/{model_prefix}_beta_tab2.rds")
      )
      tab2
    }
)

walk2(
  samples_tost, model_features$model_prefix,
  function(tost, model_prefix) {
  p1 <- TOST_figure(tost$TOST_bestpars)
  save_plot(
    filename = glue("figures/{model_prefix}_beta_tost.pdf"), p1
  )
})


walk2(
  best_si, model_features$model_prefix,
  function(si, model_prefix) {
  psi <- SI_figure(si[[2]], cowling_data)
  save_plot(
    filename = glue("figures/{model_prefix}_beta_si.pdf"), psi
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
          data = cowling_data
        )
        save_plot(
          filename = glue("figures/{model_prefix}_beta_si.pdf"), psi2
        )
    }
  }
)

x <- as.list(model_features)
x <- append(x, values = list(fits = fits))

dic <- pmap_dbl(
  x, function(mixture, recall, right_bias, model_prefix, fits) {
    samples <- rstan::extract(fits)
    DIC_alt(log_likel(samples, mixture, recall, "beta"))
  }
)

names(dic) <- x$model_prefix
## To be executed after getting Table 2s from
## various models
overall_table2 <- map_dfr(
  model_features$model_prefix,
  function(model_prefix) {
    out <- readRDS(glue("processed_stanfits/{model_prefix}_beta_tab2.rds"))
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
readr::write_csv(x, "beta_model_outputs.csv")
