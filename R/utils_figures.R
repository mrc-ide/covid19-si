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

outdir <- "processed_stanfits/release"
figs_dir <- "figures/release"


overall_table2 <- readRDS(glue("{outdir}/nf_overall_table2.rds"))
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
## Top 4 moldes
top4 <- best_unconditional[best_unconditional$model %in% overall_table2$model[1:4], ]
palette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00")
names(palette) <- overall_table2$model[1:4]



p <- ggplot(top4, aes(model, si, fill = model)) +
  gghalves::geom_half_violin(
    draw_quantiles = c(0.25, 0.5, 0.75)
    ) +
  scale_fill_manual(
    values = palette,
    breaks = overall_table2$model[1:4],
    labels = nice_model_name
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
  samples_tost, function(x) data.frame(tost = x[["TOST_bestpars"]]),
  .id = "model"
)
tost_bestpars$model <- factor(
  tost_bestpars$model, levels = overall_table2$model,
  ordered = TRUE
)
top4 <- tost_bestpars[tost_bestpars$model %in% overall_table2$model[1:4], ]
p <- ggplot(top4, aes(model, tost, fill = model)) +
  gghalves::geom_half_violin(
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) +
  scale_fill_manual(
    values = palette,
    breaks = overall_table2$model[1:4],
    labels = nice_model_name
  ) +
  theme_minimal() +
  ylab("TOST") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
  )

ggsave("{figs_dir}/nf_all_tost_distr.pdf", p)
