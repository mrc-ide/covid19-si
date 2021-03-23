process_fits <- readRDS('3a_mix_stress_testing_sim_processed_fits.rds')

######################### Plot 1. Just the parameters that are input
x <- process_fits[, c('sim', 'incubation', 'isolation', 'offset')]
x <- tidyr::gather(x, var, val, -sim)
x$sim <- factor(x$sim, levels = 1:16, ordered = TRUE)

ptable <- ggplot() +
  geom_text(data = x, aes(var, sim, label = val)) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.title = element_blank()) +
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
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(
          color = "black", fill = NA, size = 1
        ),
        strip.background = element_rect(color = "black", size = 1))

p <- ptable + pest + plot_layout(ncol = 2, widths = c(0.5, 1))
