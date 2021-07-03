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

outdir <- "processed_stanfits/relaxed_priors"
figs_dir <- "figures/relaxed_priors"


overall_table2 <- readRDS(glue("{outdir}/nf_overall_table2.rds"))
overall_table2 <- arrange(overall_table2, DIC)
names(best_si) <- model_features$model_prefix

best_unconditional <- map_dfr(
  grep("scenario3", model_features$model_prefix, value = TRUE),
  function(model) {
    s4 <- gsub("scenario3", "scenario4", model)
    data.frame(
      s3si = best_si[[model]][["unconditional"]],
      s4si = best_si[[s4]][["unconditional"]],
      s3model = model
    )
  }
)

best_unconditional$s3model <- factor(
  best_unconditional$s3model, levels = overall_table2$model,
  ordered = TRUE
)
## Top 4 moldes
top4 <- best_unconditional[best_unconditional$model %in% overall_table2$model[1:4], ]
palette <- c("#56B4E9", "#009E73", "#CC79A7", "#D55E00")
names(palette) <- c(
  "scenario3a", "scenario3a_mixture",
  "scenario3a_recall", "scenario3a_mixture_recall"
)



p <- ggplot(best_unconditional) +
  gghalves::geom_half_violin(
    aes(s3model, s3si, fill = s3model),
    draw_quantiles = c(0.25, 0.5, 0.75), side = "l"
    ) +
  gghalves::geom_half_violin(
    aes(s3model, s4si, fill = s3model),
    draw_quantiles = c(0.25, 0.5, 0.75), side = "r",
    alpha = 0.3
  ) +
  scale_fill_manual(
    values = palette,
    breaks = c(
      "scenario3a", "scenario3a_mixture",
      "scenario3a_recall", "scenario3a_mixture_recall"
    ),
    labels = c("BASELINE/ISOL", "BASELINE/ISOL + MIX",
               "BASELINE/ISOL + RECALL",
               "BASELINE/ISOL + MIX + RECALL"),
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

ggsave(glue("{figs_dir}/nf_top4_si_distr.pdf"), p)


## Similary TOST
names(samples_tost) <- model_features$model_prefix

tost_bestpars <- map_dfr(
  grep("scenario3", model_features$model_prefix, value = TRUE),
  function(model) {
    s4 <- gsub("scenario3", "scenario4", model)
    data.frame(
      s3tost = samples_tost[[model]][["TOST_bestpars"]],
      s4tost = samples_tost[[s4]][["TOST_bestpars"]],
      s3model = model
    )
  }
)

p <- ggplot(tost_bestpars) +
  gghalves::geom_half_violin(
    aes(s3model, s3tost, fill = s3model),
    draw_quantiles = c(0.25, 0.5, 0.75), side = "l"
    ) +
  gghalves::geom_half_violin(
    aes(s3model, s4tost, fill = s3model),
    draw_quantiles = c(0.25, 0.5, 0.75),
    side = "r", alpha = 0.3
  ) +
  scale_fill_manual(
    values = palette,
    breaks = c(
      "scenario3a", "scenario3a_mixture",
      "scenario3a_recall", "scenario3a_mixture_recall"
    ),
    labels = c("BASELINE/ISOL", "BASELINE/ISOL + MIX",
               "BASELINE/ISOL + RECALL",
               "BASELINE/ISOL + MIX + RECALL"),
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

ggsave(
  glue("{figs_dir}/nf_all_tost_distr.pdf"), p
)
