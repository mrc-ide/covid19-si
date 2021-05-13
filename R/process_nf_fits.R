model_features <- list(
  "mixture" = c(TRUE, FALSE),
  "left_bias" = c(TRUE, FALSE),
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

model_features$model_prefix <-ifelse(
  model_features$`left_bias`,
  glue::glue("{model_features$model_prefix}_leftbias"),
  model_features$model_prefix
)

model_features$model_prefix <-ifelse(
  model_features$`recall`,
  glue::glue("{model_features$model_prefix}_recall"),
  model_features$model_prefix
)

## Need to rename fit files, at the moment only 1 has the right
## name
index <- 16
model_features <- model_features[index, ]

if (! dir.exists("processed_stanfits")) dir.create("processed_stanfits")

fits <- map(
  model_features$model_prefix,
  function(model_prefix) {
    readRDS(glue("stanfits/{model_prefix}_nf_fit.rds"))
  }
)

table1 <- map(
  fits,
  function(fit) {
    tab1 <- fitted_params(fit)
    saveRDS(
      tab1, glue("processed_stanfits/{model_prefix}_nf_tab1.rds")
    )
    tab1
  }
)


samples_tost <- map2(
  table1, fits, function(tab1, fit) {
    estimated_TOST_nf(
      tab1, taus = seq(-20, 40, 0.1), n = 1e4,
      rstan::extract(fit)
    )
  }
)

best_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1
  ), function(tost, mixture, recall, isol, tab1) {
    estimated_SI(
      cowling_data, tost$TOST_bestpars, mixture = mixture,
      recall = recall, isol = right_bias, tab1 = tab1
    )
  }
)

mean_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1
  ), function(tost, mixture, recall, isol, tab1) {
    estimated_SI(
      cowling_data, tost$TOST_meanpars, mixture = mixture,
      recall = recall, isol = right_bias, tab1 = tab1
    )
  }
)

post_si <- pmap(
  list(
    tost = samples_tost, mixture = model_features$mixture,
    recall = model_features$recall, isol = model_features$right_bias,
    tab1 = table1
  ), function(tost, mixture, recall, isol, tab1) {
    x <- apply(
      tost$TOST_post, 2, function(inf_samples) {
        estimated_SI(
          cowling_data, inf_samples, mixture = mixture,
          recall = recall, isol = right_bias, tab1 = tab1
        )
      }
    )
    ## Extract conditional SI and make a matrix/data.frame
    map(x, ~ .[[2]]) %>% do.call(what = 'rbind')
  }
)

table2 <- pmap(
  list(
    mean_si = mean_si, best_si = best_si, post_si = post_si,
    tost = samples_tost), function(mean_si, best_si, post_si, tost) {
      samples_si <- list(
        SI_meanpars = list(SI = mean_si[[2]]),
        SI_bestpars = list(SI = best_si[[2]]),
        SI_post = post_si
      )
      ## table 2 - summary stats for sampled distributions
      tab2 <- tost_si_summary(tost, samples_si)
      saveRDS(
        tab2, glue("processed_stanfits/{model_prefix}_nf_tab2.rds")
      )
      tab2
    }
)

walk(samples_tost, function(tost) {
  p1 <- TOST_figure(tost$TOST_bestpars)
  save_plot(
    filename = glue("figures/{model_prefix}_nf_tost.png"), p1
  )
})


walk(best_si, function(si) {
  psi <- SI_figure(si[[2]], cowling_data)
  save_plot(
    filename = glue("figures/{model_prefix}_nf_si.png"), psi
  )
  psi2 <- plot_both_SIs(
    SI1 = best_si[[1]], SI2 = best_si[[2]], data = cowling_data
  )
  save_plot(
    filename = glue("figures/{model_prefix}_nf_si.png"), psi2
  )
})


