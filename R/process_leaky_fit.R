library(ggplot2)
source("R/utils.R")
source("R/utils_process_fits_common.R")
source("R/utils_process_nf_fits.R")
source("R/utils_model_selection.R")
source("R/utils_process_beta_fits.R")

fit <- readRDS("stanfits/scenario4_leaky_nf.rds")[[1]]
tab1 <- fitted_params(fit)
tost <- estimated_TOST_nf(
  tab1, taus = seq(-20, 21, 0.1), n = 1e4,
  rstan::extract(fit)
)
best_si <-  estimated_SI(
  cowling_data, tost$TOST_bestpars, mixture = FALSE,
  recall = FALSE, isol = FALSE, leaky = TRUE, tab1 = tab1
)

un_si <- unconditional_si_all(cowling_data, tost$TOST_bestpars)
pleak_par <- tab1["pleak", "best"]
leaky <- leaky_si_all(un_si, pleak_par, 1e4)

out <- map_dfr(leaky, function(x) {
  imap_dfr(
    x, function(a, b) data.frame(val = a),
    .id = "param"
  )
}, .id = "nu")

out$nu <- factor(
  out$nu, levels = as.character(as.numeric(unique(out$nu))),
  ordered = TRUE
)
x <- out[out$param %in% c("leaky_si", "not_leaky_si" ), ]
x <- x[x$nu %in% c("-4", "2", "6"), ]

ggplot() +
  ggridges::geom_density_ridges(
    data = x,
    aes(y = nu, x = val, fill = param), alpha = 0.3, col = NA
    ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())



ggplot(NULL) +
  geom_histogram(aes(cowling_data$si, y = ..density..), alpha = 0.2) +
  geom_density(aes(best_si[["unconditional"]], fill = "blue"), col = NA, alpha = 0.3) +
  geom_density(aes(best_si[["conditional"]], fill = "red"), col = NA, alpha = 0.3) +
  scale_fill_identity(
    breaks = c("blue", "red"),
    labels = c("unconditional", "conditional"),
    guide = "legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

mean_si <- estimated_SI(
  cowling_data, tost$TOST_meanpars, mixture = FALSE,
  recall = FALSE, isol = FALSE, leaky = TRUE, tab1 = tab1
)
x <- apply(
  tost$TOST_post, 2, function(inf_samples) {
    estimated_SI(
      cowling_data, inf_samples, mixture = FALSE,
      recall = FALSE, isol = FALSE, leaky = TRUE, tab1 = tab1
    )
  }
)
## Extract conditional SI and make a matrix/data.frame

post_si <- map(x, ~ .[[2]]) %>% do.call(what = 'rbind')
samples_si <- list(
  SI_meanpars = list(SI = mean_si[[2]]),
  SI_bestpars = list(SI = best_si[[2]]),
  SI_post = post_si
)
## table 2 - summary stats for sampled distributions
tab2 <- tost_si_summary(tost, samples_si)

psi <- SI_figure(best_si[[2]], cowling_data)

save_plot(
  filename = glue("figures/scenario3a_leaky_nf_si.pdf"), psi
)
