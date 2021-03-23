prefix <- "3a_mix_stress_testing_sim"
process_fits <- readRDS('stanfits/{prefix}_processed_fits.rds')

######################### Plot 1. Just the parameters that are input
x <- process_fits[, c('sim', 'incubation', 'isolation', 'offset')]
x <- tidyr::gather(x, var, val, -sim)
x$sim <- factor(x$sim, levels = 1:16, ordered = TRUE)

ptable <- ggplot(x) +
  geom_text(aes(1, sim, label = val)) +
  facet_wrap(~var, ncol = 3) +
  theme_minimal() +
  theme(
    axis.title = element_blank(), axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  theme(
    panel.spacing = unit(.05, "lines"),
    panel.border = element_rect(
      color = "black", fill = NA, size = 1
    ),
    strip.background = element_rect(color = "black", size = 1))


######################### Plot 2. Estimated parameters
y <- select(
  process_fits, sim, `mu_2.5%`:`mu_97.5%`, `sd_2.5%`:`sd_97.5%`,
  `pinvalid_2.5%`:`pinvalid_97.5%`,
  mu_true = true_mean, sd_true = true_sd, pinvalid_true = pinvalid
  )
####### Reorganise, should have set-up like this in the first place!!
y <- tidyr::gather(y, var, val, -sim)
y <- tidyr::separate(y, col = "var", into = c("var", "qntl"), sep = "_")
y$qntl[is.na(y$qntl)] <- "true_val"
y <- tidyr::spread(y, qntl, val)
y$sim <- factor(y$sim, levels = 1:16, ordered = TRUE)

pest <- ggplot(y) +
  geom_point(aes(`50%`,sim)) +
  geom_linerange(aes(xmin = `25%`, xmax = `75%`, y = sim)) +
  geom_point(aes(`true`, sim), shape = 4) +
  facet_wrap(~var, ncol = 3, scales = 'free_x') +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(axis.title = element_blank()) +
  theme(panel.spacing = unit(.1, "lines"),
        panel.border = element_rect(
          color = "black", fill = NA, size = 1
        ),
        strip.background = element_rect(color = "black", size = 1))

p <- ptable + pest + plot_layout(ncol = 2, widths = c(0.5, 1))

cowplot::save_plot(glue::glue("figures/{prefix}all_vars.pdf"), p)

outfiles <- glue::glue("data/{prefix}_{seq_along(mixed)}data.rds")
training <- map(outfiles, readRDS)

outfiles <- glue::glue("stanfits/{prefix}_{seq_along(mixed)}posterior.rds")
posterior <- map(outfiles, readRDS)

plots <- map2(
  training, posterior, function(trng, post) {
    ggplot() +
      geom_density(
        aes(trng$si, fill = 'red'),
        col = NA, alpha = 0.3
      ) +
      geom_density(
        aes(post$si, fill = 'blue'),
        col = NA, alpha = 0.3
      ) +
      scale_fill_identity(
        breaks = c('red', 'blue'), labels = c('Training', 'Posterior'),
        guide = 'legend'
      ) +
      theme_minimal() +
      theme(
        legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank()
      )
  }
)



p <- wrap_plots(
  plots, ncol = 4, nrow = 4, guides = 'collect', byrow = TRUE
) &
  theme(legend.position = 'top')
## cowplot defaults don't look good in this case.
ggsave(glue::glue("figures/{prefix}all_si.pdf"), p)
