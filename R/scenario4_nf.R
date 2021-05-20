prefix <- "scenario4a_nf"
infile <- glue::glue("stan-models/{prefix}.stan")
params_inc <- params_real$inc_par2
si_vec <- seq(-20, max_valid_si)

fit_4a <- stan(
  file = here::here(infile),
  data = list(
    N = nrow(cowling_data),
    si = cowling_data$si,
    nu = cowling_data$nu,
    max_shed = max_shed,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    M = length(si_vec),
    si_vec = si_vec,
    width = 0.1
  ),
  verbose = TRUE
)

saveRDS(
  fit_4a, glue::glue("stanfits/{prefix}_fit.rds")
)

ggmcmc(
  ggs(fit_4a), glue::glue("checks/{prefix}_check.pdf")
)

