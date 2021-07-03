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

outdir <- "processed_stanfits/discrete_pairs"
figs_dir <- "figures/discrete_pairs"


overall_table2 <- readRDS(glue("{outdir}/nf_overall_table2.rds"))
overall_table2 <- arrange(overall_table2, DIC)
names(best_si) <- model_features$model_prefix

best_unconditional <- map_dfr(
  overall_table2$model[1:4],
  function(model) {
    s3 <- ifelse(
      grepl('scenario3', model),
      model,
      gsub('scenario4', 'scenario3', model)
    )
    s4 <- ifelse(
      grepl('scenario3', model),
      gsub('scenario4', 'scenario3', model),
      model
    )
    data.frame(
      s3si = best_si[[s3]][["unconditional"]],
      s4si = best_si[[s4]][["unconditional"]],
      model = model
    )
  }
)

best_unconditional$model <- factor(
  best_unconditional$model, levels = overall_table2$model[1:4],
  ordered = TRUE
)
## Top 4 models
##s3models <- grep('scenario3', overall_table2$model, value = TRUE)
##top4 <- best_unconditional[best_unconditional$model %in% overall_table2$model[1:4], ]
palette <- c("#56B4E9", "#009E73", "#CC79A7", "#D55E00")
names(palette) <- overall_table2$model[1:4]


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
    breaks = overall_table2$model[1:4],
    labels = nice_model_name(overall_table2$model[1:4]),
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
  overall_table2$model[1:4],
  function(model) {
    s3 <- ifelse(
      grepl('scenario3', model),
      model,
      gsub('scenario4', 'scenario3', model)
    )
    s4 <- ifelse(
      grepl('scenario3', model),
      gsub('scenario4', 'scenario3', model),
      model
    )
    data.frame(
      s3tost = samples_tost[[s3]][["TOST_bestpars"]],
      s4tost = samples_tost[[s4]][["TOST_bestpars"]],
      model = model
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
    breaks = overall_table2$model[1:4],
    labels = nice_model_name(overall_table2$model[1:4]),
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
