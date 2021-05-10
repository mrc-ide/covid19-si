## Scenario 3, no mixture
prefix <- "scenario3a_nf"
infile <- glue::glue("stan-models/{prefix}.stan")
params_inc <- params_real$inc_par2
si_vec <- seq(-20, max_valid_si)

fit_3a <- stan(
  file = here::here(infile),
  data = list(
    N = nrow(cowling_data),
    si = cowling_data$si,
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
  fit_3a, glue::glue("stanfits/{prefix}_fit.rds")
)

ggmcmc(
  ggs(fit_3a), glue::glue("checks/{prefix}_check.pdf")
)
