param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2", "inf_par3"),
  params_inc = c("inc_par1", "inc_par2"),
  params_pinv = c("pinvalid1", "pinvalid2", "pinvalid3"),
  params_recall = c("recall1", "recall2", "recall3"),
  stringsAsFactors = FALSE
)

## Pick a large number here since after filtering, you will
## be left with a much smaller data set.

nsim <- 10000
invalid_si <- rbeta(
  nsim, shape1 = alpha_invalid, shape2 = beta_invalid
)

invalid_iso <- rgamma(
  nsim, shape = params_iso$shape, scale = params_iso$scale
)
invalid_si <- data.frame(si = invalid_si, nu = invalid_iso)


valid_unfiltered <- pmap(
  param_grid,
  function(params_inf, params_inc, params_pinv, params_recall) {
    outfile <- paste(
      params_inf, params_inc, params_pinv, sep = "_"
    )

    params_inf <- params[[params_inf]]
    params_inc <- params[[params_inc]]

    params_inf <- beta_muvar2shape1shape2(
      params_inf$mean_inf/max_shed, params_inf$sd_inf^2 / max_shed^2
    )
    valid <- simulate_si(
      params_inc$mean_inc, params_inc$sd_inc, params_inf$shape1,
      params_inf$shape2, max_shed, mean_iso, sd_iso, nsim
    )
    valid
  }
)

valid_filtered <- map(
  valid_unfiltered, function(simulated_si) {
    simulated_si[simulated_si$t_1 <= simulated_si$nu, ]
 }
)

invalid_remapped <- map(
  valid_unfiltered,
  function(valid) {
    df <- invalid_si
    max_si <- max(valid$si)
    df$si <- (max_si - min_si) * df$si + min_si
    df
  }
)

mixed <- pmap(
  list(
    valid = valid_filtered,
    invalid = invalid_remapped,
    params_pinv = param_grid$params_pinv
  ),
  function(valid, invalid, params_pinv) {
    pinvalid <- params[[params_pinv]]
    toss <- runif(nrow(valid))
    rbind(
      valid[toss > pinvalid , c("si", "nu")],
      invalid[toss <= pinvalid ,c("si", "nu")]
    )
  }
)


with_recall_bias <- map2(
  mixed,
  param_grid$params_recall,
  function(df, recall_true) {
    recall_true <- params[[recall_true]]
    df$p_si <- exp(abs(df$si - df$nu) * -recall_true)
    idx <- sample(
      nrow(df), nrow(df), replace = TRUE, prob = df$p_si
    )
    df[idx, ]
  }
)

all_data <- pmap(
  list(
    x = valid_unfiltered,
    y = valid_filtered,
    a = mixed,
    b = with_recall_bias
  ),
  function(x, y, a, b) {
    rbind(
      data.frame(type = "all valid", si = x$si, nu = x$nu),
      data.frame(type = "valid filtered", si = y$si, nu = y$nu),
      data.frame(type = "mixed", si = a$si, nu = a$nu),
      data.frame(type = "mixed and recall", si = b$si, nu = b$nu)
    )
  }
)

iwalk(
  all_data,
  function(df, i) saveRDS(df, glue::glue("data/simulated_{i}.rds"))
)

palette <- c(
 "all valid" = "#999999",
 "valid filtered" = "#E69F00",
  "mixed" = "#56B4E9",
  "mixed and recall" = "#CC79A7"
)

iwalk(
  all_data,
  function(df, index) {
    p <- ggplot(df, aes(si, fill = type)) +
      geom_density(alpha = 0.4, col = NA) +
      scale_fill_manual(
        values = palette,
        breaks = c("all valid", "valid filtered",
                   "mixed", "mixed and recall"),
        labels = c("All possible SIs",
                   "Infection of infectee before isolation of primary",
                   "Mixture distribution",
                   "Mixture distribution sampled with recall bias"),
        guide = guide_legend(nrow = 2)
      ) +
      theme_minimal() +
      theme(legend.position = "top", legend.title = element_blank()) +
      xlab("Serial Interval") +
      ylab("Density")

    ggsave(glue::glue("figures/{index}.pdf"), p)
  }
)

iwalk(
  all_data,
  function(df, index) {
    p <- ggplot(df, aes(nu, si, col = type)) +
      geom_point() +
      scale_color_manual(
        values = palette,
        breaks = c("all valid", "valid filtered",
                   "mixed", "mixed and recall"),
        labels = c("All possible SIs",
                   "Infection of infectee before isolation of primary",
                   "Mixture distribution",
                   "Mixture distribution sampled with recall bias"),
        guide = guide_legend(nrow = 2)
      ) +
      theme_minimal() +
      theme(legend.position = "top", legend.title = element_blank()) +
      xlab("Delay from symptom onset to isolation") +
      ylab("Serial Interval")

    ggsave(glue::glue("figures/nu_vs_si_{index}.pdf"), p)
  }
)





