## Compare conditional distributions across
## TOST
compare_across_distrs <- function(obs, nf, sn, gm) {
  alpha <- 0.2
  p <- ggplot() +
    geom_histogram(
      data = obs, aes(si, y = ..density.., fill = "gray77"),
      alpha = 0.8,
      binwidth = 1
    ) +
    geom_density(
      aes(nf, fill = "red"), alpha = alpha, colour = NA
    ) +
    geom_density(
      aes(sn, fill = "blue"), alpha = alpha, colour = NA
    ) +
    geom_density(
      aes(gm, fill = "green"), alpha = alpha, colour = NA
    ) +
    scale_fill_identity(
      guide = "legend",
      labels = c("Data", "NF", "Skew normal", "gamma"),
      breaks = c("gray77", "red", "blue", "green")
    ) +
    theme_minimal() +
    xlab("Serial Interval (days)") +
    theme(legend.title = element_blank()) +
    theme(legend.position = "bottom")
  p
}

obs_data <- data_discrete_pairs_s3_s4mix
meta <- "s3s4pairs"
##nf_dir <- "processed_stanfits/maxshed21_nfpriors/discrete_pairs"
gamma_dir <- "processed_stanfits/gamma/s3s4pairs"
sn_dir <- "processed_stanfits/skew_normal/s3s4pairs"

nf_si <- readRDS(glue("{nf_dir}/best_si_nf.rds"))
sn_si <- readRDS(glue("{sn_dir}/best_si_skew_normal.rds"))
gamma_si <- readRDS(glue("{gamma_dir}/best_si_gamma.rds"))

outdir <- glue("figures/compare_all/{meta}")
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
pwalk(
  list(
    nf = nf_si, sn = sn_si, gm = gamma_si,
    prefix = model_features$model_prefix[model_features$right_bias]
  ),
  function(nf, sn, gm, prefix) {
    x <- nf[["conditional"]]
    y <- sn[["conditional"]]
    z <- gm[["conditional"]]
    p <- compare_across_distrs(
      obs_data, x, y, x
    )
    p <- p + ggtitle(snakecase::to_title_case(prefix))
    ggsave(glue("{outdir}/{prefix}_all_distrs.png"), p)
  }
)



## Table 2.
outdir <- glue("processed_stanfits/compare_all/{meta}")
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
nf_t2 <- readRDS(glue("{nf_dir}/nf_overall_table2.rds"))
nf_t2$TOST <- "NF"
sn_t2 <- readRDS(glue("{sn_dir}/skew_normal_overall_table2.rds"))
sn_t2$TOST <- "skew normal"
gamma_t2 <- readRDS(glue("{gamma_dir}/gamma_overall_table2.rds"))
gamma_t2$TOST <- "gamma"


overall <- rbind(sn_t2, gamma_t2)
overall$model <- nice_model_name(overall$model)
##overall$model <- paste("S3S4", overall$model)
saveRDS(overall, glue("{outdir}/table2_all_distrs.rds"))
