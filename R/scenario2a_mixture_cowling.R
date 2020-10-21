data <- readRDS("data/cowling_data_clean.rds")
data_pos <- data%>%
  filter(si > 0)%>%
  filter(onset_first_iso > 0)%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)


fit_mixture <- stan(
  file = here::here("stan-models/scenario2a_mixture.stan"),
  data = list(
    N = length(data_pos$si),
    si = data_pos$si,
    nu = data_pos$nu,
    max_shed = max_shed,
    alpha2 = params_inc_og[["shape"]],
    beta2 = 1 / params_inc_og[["scale"]],
    alpha_invalid = alpha_invalid,
    beta_invalid = beta_invalid,
    max_si = max(data_pos$si) + 0.001,
    min_si = min(data_pos$si),
    width = min(data_pos$si) / 2
  ),
  chains = 2,
  iter = 5000,
  verbose = TRUE
  ## control = list(adapt_delta = 0.99)
)
