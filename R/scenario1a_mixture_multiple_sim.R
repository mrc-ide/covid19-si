## Set up params grid
param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2", "inf_par3"),
  params_inc = c("inc_par1", "inc_par2"),
  params_pinv = c("pinvalid1", "pinvalid2", "pinvalid3"),
  stringsAsFactors = FALSE
)

nsim <- 500
alpha_invalid <- 0.5
beta_invalid <- 0.5
min_si <- -2

out <- pmap(
  param_grid,
  function(params_inf, params_inc, params_pinv) {
    outfile <- paste(
      params_inf, params_inc, params_pinv, sep = "_"
    )

    params_inf <- params[[params_inf]]
    params_inc <- params[[params_inc]]
    pinvalid <- params[[params_pinv]]
    params_inf <- beta_muvar2shape1shape2(
      params_inf$mean_inf/max_shed, params_inf$sd_inf^2 / max_shed^2
    )
    simulated_si <- simulate_1a_mix(
      params_inc$mean_inc, params_inc$sd_inc, params_inf, pinvalid, nsim
    )

    p <- plot_1a_mix(
      simulated_si[[1]], simulated_si[[2]], simulated_si[[3]]
    )
    ggsave(glue::glue("figures/{outfile}.png"), p)

    simulated_si <- simulated_si[[3]]

    saveRDS(simulated_si, glue::glue("stanfits/simulated_{outfile}.rds"))

    width <- min(simulated_si[simulated_si > 0]) / 2
    params2_inc <- gamma_mucv2shapescale(
      params_inc$mean_inc, params_inc$sd_inc/ params_inc$mean_inc
    )
    fit_mixture <- stan(
      file = here::here("stan-models/scenario1a_mixture_general.stan"),
      data = list(
        N = length(simulated_si),
        si = simulated_si,
        max_shed = max_shed,
        alpha2 = params2_inc[["shape"]],
        beta2 = 1 / params2_inc[["scale"]],
        alpha_invalid = alpha_invalid,
        beta_invalid = beta_invalid,
        max_si = max(simulated_si) + 0.001,
        min_si = min(simulated_si) - 0.001,
        width = width
      ),
      chains = 5,
      iter = 2000,
      verbose = TRUE
      ## control = list(adapt_delta = 0.99)
    )
    saveRDS(fit_mixture, glue::glue("stanfits/{outfile}.rds"))

    best_params <- map_estimates(fit_mixture)
    posterior_si <- simulate_1a_mix(
      params_inc$mean_inc, params_inc$sd_inc,
      list(shape1 = best_params[["alpha1"]],
           shape2 = best_params[["beta1"]]),
      best_params[["pinvalid"]], 10000
    )
    saveRDS(
      posterior_si, glue::glue("stanfits/posterior_si_{outfile}.rds")
    )

    p2 <- ggplot() +
      geom_density(aes(simulated_si, fill = "blue"), alpha = 0.3) +
      geom_density(aes(posterior_si[[3]], fill = "red"), alpha = 0.3) +
      scale_fill_identity(
        guide = "legend",
        labels = c("Simulated", "Posterior"),
        breaks = c("blue", "red")
      ) +
      xlab("Serial Interval") +
      ylab("Probability Density") +
      theme_minimal() +
      theme(legend.title = element_blank())

    ggsave(glue::glue("figures/posterior_{outfile}.png"), p2)

  }
)

pinvalid_compare <- pmap(
  param_grid,
  function(params_inf, params_inc, params_pinv) {
    infile <- paste(
      params_inf, params_inc, params_pinv, sep = "_"
    )
    fit <- readRDS(glue::glue("stanfits/{infile}.rds"))
    fit_params <-  rstan::extract(fit)
    out <- map_dfr(fit_params, quantile_as_df, .id = "param")
    saveRDS(out, glue::glue("stanfits/fitted_params_{infile}.rds"))
   }
)

compare <- pmap(
  param_grid,
  function(params_inf, params_inc, params_pinv) {
    infile <- paste(
      params_inf, params_inc, params_pinv, sep = "_"
    )
    simulated <- readRDS(
      glue::glue("stanfits/simulated_{infile}.rds")
    )
    posterior <- readRDS(
      glue::glue("stanfits/posterior_si_{infile}.rds")
    )
    simulated_summary <- quantile_as_df(simulated)
    simulated_summary$category <- "Simulated"
    posterior_summary <- quantile_as_df(posterior[[3]])
    posterior_summary$category <- "Posterior"
    rbind(simulated_summary, posterior_summary)
})

compare <- dplyr::bind_rows(compare, .id = "sim")

compare <- tidyr::spread(compare, var, val)
compare$sim <- factor(compare$sim, levels = 1:18, ordered = TRUE)

p <- ggplot(compare) +
  geom_point(
    aes(sim, `50%`, col = category),
    size = 2,
    position = position_dodge(width = 0.3)
  ) +
  geom_linerange(
    aes(sim, ymin = `2.5%`, ymax = `97.5%`, col = category),
    size = 1.1,
    position = position_dodge(width = 0.3)
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank()) +
  xlab("Simulation") +
  ylab("Serial Interval (Median and 95% CrI)")

ggsave("figures/scenario1a_mix_multiple_sim.png", p)
