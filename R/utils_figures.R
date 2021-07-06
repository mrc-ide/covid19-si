nice_model_name <- function(model) {
  prefix <- case_when(
    grepl("scenario3", model) ~ "BASELINE",
    grepl("scenario4", model) ~ "ISOL"
  )
  suffix1 <- case_when(
    grepl("mixture", model) ~ " + MIX",
    TRUE ~ ""
  )
  suffix2 <- case_when(
    grepl("recall", model) ~ " + RECALL",
    TRUE ~ ""
  )

  paste0(prefix, suffix1, suffix2)
}

palette <- c(
  scenario3a = "#56B4E9", scenario4a = "#56B4E9",
  scenario3a_mixture = "#009E73",
  scenario4a_mixture = "#009E73",
  scenario3a_recall = "#CC79A7",
  scenario4a_recall = "#CC79A7",
  scenario3a_mixture_recall = "#D55E00",
  scenario4a_mixture_recall = "#D55E00"
)


meta_model <- "maxshed21_nfpriors/s3s4mix"
fit_dir <- glue("stanfits/{meta_model}")
outdir <- glue("processed_stanfits/{meta_model}")
figs_dir <- glue("figures/{meta_model}")

overall_table2 <- readRDS(glue("{outdir}/nf_overall_table2.rds"))
overall_table2 <- arrange(overall_table2, DIC)

best_si <- readRDS(glue("{outdir}/best_si_nf.rds"))
names(best_si) <- model_features$model_prefix

s3models <- c(
  "scenario3a", "scenario3a_mixture",
  "scenario3a_recall", "scenario3a_mixture_recall"
)

best_unconditional <- map_dfr(
  s3models,
  function(s3) {
    s4 <- gsub('scenario3', 'scenario4', s3)
    data.frame(
      s3si = best_si[[s3]][["unconditional"]],
      s4si = best_si[[s4]][["unconditional"]],
      model = s3
    )
  }
)

best_unconditional$model <- factor(
  best_unconditional$model,
  levels = s3models,
  ordered = TRUE
)


p <- ggplot(best_unconditional) +
  gghalves::geom_half_violin(
    aes(model, s3si, fill = model),
    draw_quantiles = c(0.25, 0.5, 0.75), side = "l"
  ) +
  gghalves::geom_half_violin(
    aes(model, s4si, fill = model),
    draw_quantiles = c(0.25, 0.5, 0.75), side = "r",
    alpha = 0.3
  ) +
  scale_fill_manual(
    values = palette,
    breaks = levels(best_unconditional$model),
    labels = nice_model_name,
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_minimal() +
  ylab("Serial Interval") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
  )

ggsave(glue("{figs_dir}/nf_top4_si_distr.png"), p)


## Similary TOST
samples_tost <- readRDS(glue("{outdir}/samples_tost_nf.rds"))
names(samples_tost) <- model_features$model_prefix

tost_bestpars <- map_dfr(
  s3models,
  function(s3) {
    s4 <- gsub('scenario3', 'scenario4', s3)
    data.frame(
      s3tost = samples_tost[[s3]][["TOST_bestpars"]],
      s4tost = samples_tost[[s4]][["TOST_bestpars"]],
      model = s3
    )
  }
)

p <- ggplot(tost_bestpars) +
  gghalves::geom_half_violin(
    aes(model, s3tost, fill = model),
    draw_quantiles = c(0.25, 0.5, 0.75), side = "l"
    ) +
  gghalves::geom_half_violin(
    aes(model, s4tost, fill = model),
    draw_quantiles = c(0.25, 0.5, 0.75),
    side = "r", alpha = 0.3
  ) +
  scale_fill_manual(
    values = palette,
    breaks = levels(best_unconditional$model),
    labels = nice_model_name,
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_minimal() +
  ylab("TOST") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
  )

ggsave(glue("{figs_dir}/nf_all_tost_distr.png"), p)
