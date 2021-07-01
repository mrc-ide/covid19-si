fit_dir <- "stanfits/release"
outdir <- "processed_stanfits/release"
figs_dir <- "figures/release"

index <- map_lgl(
  model_features$model_prefix,
  function(model_prefix) file.exists(glue("{fit_dir}/{model_prefix}_nf_fit.rds"))
)

model_features <- model_features[index, ]

if (! dir.exists("processed_stanfits")) dir.create("processed_stanfits")
check <- "\U2713"
fits <- map(
  model_features$model_prefix,
  function(model_prefix) {
    message("Reading fit file for ", model_prefix)
    readRDS(glue("{fit_dir}/{model_prefix}_nf_fit.rds"))
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
      tab1, glue("{outdir}/{model_prefix}_nf_tab1.rds")
    )
    tab1
  }
)

## table1 <- map(
##   model_features$model_prefix,
##   function(x) readRDS(glue("{outdir}/{x}_nf_tab1.rds"))
## )

samples_tost <- map2(
  table1, fits, function(tab1, fit) {
    estimated_TOST_nf(
      tab1, taus = seq(-20, 40, 0.1), n = 1e4,
      rstan::extract(fit)
    )
  }
)

saveRDS(samples_tost, "{outdir}/samples_tost_nf.rds")

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
      recall = recall, isol = isol, leaky = FALSE, tab1 = tab1
    )
  }
)

saveRDS(best_si, "{outdir}/best_si_nf.rds")

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
      tost$TOST_post, 2, function(inf_samples) {
        estimated_SI(
          cowling_data, inf_samples, mixture = mixture,
          recall = recall, isol = isol, tab1 = tab1
        )
      }
    )
    ## Extract conditional SI and make a matrix/data.frame
    map(x, ~ .[["unconditional"]]) %>% do.call(what = 'rbind')
  }
)

saveRDS(post_si, "{outdir}/nf_post_si.rds")
## post_si <- readRDS("{outdir}/nf_post_si.rds")

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
        tab2, glue("{outdir}/{model_prefix}_nf_tab2.rds")
      )
      tab2
    }
)

walk2(
  samples_tost, model_features$model_prefix,
  function(tost, model_prefix) {
  p1 <- TOST_figure(tost$TOST_bestpars)
  ggsave(
    filename = glue("{figs_dir}/{model_prefix}_nf_tost.pdf"), p1
  )
})


walk2(
  best_si, model_features$model_prefix,
  function(si, model_prefix) {
  psi <- SI_figure(si[[2]], cowling_data)
  ggsave(
    filename = glue("{figs_dir}/{model_prefix}_nf_si.pdf"), psi
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
        ggsave(
          filename = glue("{figs_dir}/{model_prefix}_nf_si.pdf"), psi2
        )
    }
  }
)

x <- as.list(model_features)
x <- append(x, values = list(fits = fits))

dic <- pmap_dbl(
  x, function(mixture, recall, right_bias, model_prefix, fits) {
    samples <- rstan::extract(fits)
    DIC_alt(log_likel(samples, mixture, recall, "nf"))
  }
)
names(dic) <- x$model_prefix
## To be executed after getting Table 2s from
## various models
overall_table2 <- map_dfr(
  model_features$model_prefix,
  function(model_prefix) {
    out <- readRDS(glue("{outdir}/{model_prefix}_nf_tab2.rds"))
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
saveRDS(x, "{outdir}/nf_overall_table2.rds")



## For manuscript
for_ms <- select(x, model, mean_inf, sd_inf, mean_si, sd_si, DIC)
for_ms <- arrange(for_ms, DIC)
for_ms <- left_join(
  for_ms, model_features, by = c("model" = "model_prefix")
)
for_ms <- select(
 for_ms, mixture, recall, right_bias, mean_inf, sd_inf, mean_si, sd_si, DIC
)
## stargazer::stargazer(for_ms, summary = FALSE, rownames = FALSE)
overall_table2 <- readRDS("{outdir}/nf_overall_table2.rds")
overall_table2 <- arrange(overall_table2, DIC)
names(best_si) <- model_features$model_prefix
## Reorder best_si
best_unconditional <- map_dfr(
  best_si, function(x) {
    data.frame(si = x[["unconditional"]])
  }, .id = "model"
)

best_unconditional$model <- factor(
  best_unconditional$model, levels = overall_table2$model,
  ordered = TRUE
)

p <- ggplot(best_unconditional, aes(model, si)) +
  gghalves::geom_half_violin(
    fill = "red", alpha = 0.3,
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) +
  ##geom_boxplot(width = 0.1, color = "grey", alpha = 0.5, coef = 1) +
  theme_minimal() +
  ylab("Serial Interval") +
  scale_x_discrete(
    breaks = c("scenario4a_mixture", "scenario3a_mixture_recall", "scenario4a_mixture_recall",
               "scenario4a_recall", "scenario3a_mixture", "scenario3a_recall"),
    labels = c("ISOL + MIX", "BASELINE + MIX + RECALL", "ISOL + MIX + RECALL",
               "ISOL + RECALL", "BASELINE + MIX", "BASELINE + RECALL")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )

ggsave("{figs_dir}/nf_all_si_distr.pdf", p)


## Similary TOST
names(samples_tost) <- model_features$model_prefix
tost_bestpars <- map_dfr(
  samples_tost, function(x) data.frame(tost = x[["TOST_bestpars"]]),
  .id = "model"
)
tost_bestpars$model <- factor(
  tost_bestpars$model, levels = overall_table2$model,
  ordered = TRUE
)

p <- ggplot(tost_bestpars, aes(model, tost)) +
  gghalves::geom_half_violin(
    fill = "red", alpha = 0.3,
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) +
  ##geom_boxplot(width = 0.1, color = "grey", alpha = 0.5, coef = 1) +
  theme_minimal() +
  ylab("TOST") +
  scale_x_discrete(
    breaks = c("scenario4a_mixture", "scenario3a_mixture_recall", "scenario4a_mixture_recall",
               "scenario4a_recall", "scenario3a_mixture", "scenario3a_recall"),
    labels = c("ISOL + MIX", "BASELINE + MIX + RECALL", "ISOL + MIX + RECALL",
               "ISOL + RECALL", "BASELINE + MIX", "BASELINE + RECALL")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )

ggsave("{figs_dir}/nf_all_tost_distr.pdf", p)
