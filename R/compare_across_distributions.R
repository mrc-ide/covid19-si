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

obs_data <- data_s3_s4mix
meta <- "s3s4"
nf_dir <- "processed_stanfits/maxshed21_nfpriors/s3s4mix"
gamma_dir <- "processed_stanfits/gamma/s3s4"
sn_dir <- "processed_stanfits/skew_normal/s3s4mix"

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



