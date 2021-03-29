prefix <- "4a_mix"
process_fits <- readRDS(glue::glue('stanfits/{prefix}_processed_fits.rds'))

######################### Plot 1. Just the parameters that are input
x <- process_fits[, c('sim', 'incubation', 'isolation', 'offset')]
x <- tidyr::gather(x, var, val, -sim)
nsims <- length(unique(x$sim))
x$sim <- factor(x$sim, levels = seq_len(nsims), ordered = TRUE)

## For a grid with large number of rows, we want to split them
## over multiple pages. nrow_page is the number of rows in a single
## page.
nrow_page <- 18
## This will always be 3, for the three variables we have.
ncol_page <- 3
## Make sure nrow_page is a divisor of npages
npages <- nsims / nrow_page

ptables <- map(
  1:npages, function(page) {
   ggplot(x) +
     geom_text(aes(1, 0.5, label = val), size = 3) +
     ggforce::facet_grid_paginate(
       sim~var, nrow = nrow_page, ncol = ncol_page, page = page, switch = "y"
       ) +
     ylim(0, 1) +
     theme_minimal() +
     theme(
       axis.title = element_blank(), axis.ticks = element_blank(),
       axis.text = element_blank(),
       panel.spacing = unit(.1, "lines"),
       panel.border = element_rect(
      color = "black", fill = NA, size = 0.5
      ),
      strip.background = element_rect(color = "black", size = 1)
     )
  }
)


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
y$sim <- factor(y$sim, levels = seq_len(nsims), ordered = TRUE)

y$ylevel <- 0.5 ##as.integer(y$sim) + 0.5

pest <- map(
  seq_len(npages),
  function(page) {
    ggplot(y) +
      geom_point(aes(`50%`, ylevel)) +
      geom_linerange(aes(xmin = `25%`, xmax = `75%`, y = ylevel)) +
      geom_point(aes(`true`, ylevel), shape = 4) +
      ggforce::facet_grid_paginate(
        sim~var, nrow = 18, ncol = ncol_page, page = page,
        scales = 'free_x'
        ) +
      ylim(0, 1) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing = unit(.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1)
      )
  }
)

plots <- map2(
  ptables, pest, function(p1, p2) {
   p <- wrap_plots(p1, p2, ncol = 2, widths = c(0.5, 1))
   p
  }
)

iwalk(plots, function(p, i) {
  ggsave(glue::glue("figures/{prefix}_all_vars_{i}.png"), p)
})



outfiles <- glue::glue("data/{prefix}_{levels(x$sim)}_data.rds")
training <- map(outfiles, readRDS)

outfiles <- glue::glue("stanfits/{prefix}_{levels(x$sim)}_posterior.rds")
posterior <- map(outfiles, readRDS)

training <- dplyr::bind_rows(training, .id = "sim")
posterior <- dplyr::bind_rows(posterior, .id = "sim")
training$category <- "Training"
posterior$category <- "Posterior"
x <- rbind(training, posterior)

x$sim <- factor(x$sim, levels = seq_len(nsims), ordered = TRUE)

nrow_page <- 9
ncol_page <- 4
npages <- nsims / (nrow_page * ncol_page)

walk(
  seq_len(npages), function(page) {
    p <- ggplot(x) +
      geom_density(
        aes(si, fill = category),
        col = NA, alpha = 0.3
      ) +
      geom_density(
        aes(si, fill = category),
        col = NA, alpha = 0.3
      ) +
      ggforce::facet_wrap_paginate(~sim, nrow = nrow_page, ncol = ncol_page, page = page, scales = "free_y") +
      theme_minimal() +
      theme(
        legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank()
      )
    ggsave(glue::glue("figures/{prefix}_all_si_{page}.png"), p)
  }
)







