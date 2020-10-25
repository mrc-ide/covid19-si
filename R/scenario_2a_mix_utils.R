## Switiching to density as valid and invalid SI might have very
## different number of rows
plot_2a_mix <- function(valid_si, invalid_si, simulated_si) {

  p1 <- ggplot() +
    geom_density(aes(valid_si$si, fill = "blue"), alpha = 0.3) +
    geom_density(aes(invalid_si$si, fill = "red"), alpha = 0.3) +
    theme_minimal() +
    scale_fill_identity(
      guide = "legend",
      breaks = c("blue", "red"),
      labels = c("Valid", "Invalid")
    ) +
    xlim(0, NA) +
    xlab("Serial Interval") +
    theme(legend.title = element_blank())

  p2 <- ggplot() +
    geom_density(aes(simulated_si$si), fill = "purple", alpha = 0.3) +
    theme_minimal() +
    xlim(0, NA) +
    xlab("Serial Interval")


  p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "hv")

  p
}

simulate_2a_mix <- function(mean_inc, sd_inc, params_inf, mean_iso,
                            sd_iso, pinvalid, nsim) {

  params_iso <- epitrix::gamma_mucv2shapescale(mean_iso,
                                               sd_iso / mean_iso)
  ## Generate a large number of pairs because we need to filter
  ## on t_1 < nu
  nsim_local <- 10 * nsim
  valid_si <- simulate_si(
    mean_inc, sd_inc, params_inf$shape1, params_inf$shape2,
    max_shed, mean_iso, sd_iso, nsim = nsim_local
  )
  valid_si$nu <- rgamma(
    nsim_local, shape = params_iso$shape, scale = params_iso$scale
  )
  valid_si <- valid_si[valid_si$t_1 < valid_si$nu, ]

  max_si <- max(valid_si$si)

  invalid_si <- (max_si - min_si)*
    rbeta(nsim_local, shape1 = alpha_invalid, shape2 = beta_invalid)
  invalid_si <- invalid_si + min_si
  invalid_iso <- rgamma(
    nsim_local, shape = params_iso$shape, scale = params_iso$scale
  )
  invalid_si <- data.frame(si = invalid_si, nu = invalid_iso)


  simulated_si <- valid_si[0, ]

  for (idx in seq_len(nsim)) {
    toss <- runif(1)
    if (toss <= pinvalid) {
      simulated_si <- rbind(
        simulated_si, invalid_si[sample(nrow(invalid_si), 1), ]
      )
    } else {
      simulated_si <- rbind(
        simulated_si, valid_si[sample(nrow(valid_si), 1), c("si", "nu")]
      )
    }
  }

  list(
    valid_si = valid_si,
    invalid_si = invalid_si,
    simulated_si = simulated_si
  )
}

