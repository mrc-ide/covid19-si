nf_pdf <- function(t, a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
  numerator <- a + b
  growing <- a * exp(b * (t - tmax))
  failing <- b * exp(-a * (t - tmax))
  (numerator / (growing + failing))^c
}



# run scenario 3 on the cowling data

# offset - must be a negative number!
offset <- -20

# load the data
data <- readRDS("data/cowling_data_clean.rds")

# sub-set to only incude those SIs that are possible under our assumed offset
data_offset <- data %>%
  filter(si > offset) %>%
  filter(onset_first_iso > offset) %>%
  mutate(si = as.numeric(si)) %>%
  dplyr::rename(nu = onset_first_iso)

# fit the model
si_vec <- seq(-20, 40, 1)

fits_4a <- stan(
  file = here::here("stan-models/scenario4a_mixture.stan"),
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    nu = data_offset$nu,
    max_shed = 21,
    alpha2 = params_real$inc_par2[["shape"]],
    beta2 = 1 / params_real$inc_par2[["scale"]],
    alpha_invalid = 1,
    beta_invalid = 1,
    max_valid_si = 40,
    min_valid_si = -20,
    max_invalid_si = 40,
    min_invalid_si = -20,
    width = 0.1,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
    ##tmax = 0
  ),
  chains = 2, iter = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


out <- rstan::extract(fits_4a)
idx <- which.max(out[["lp__"]])
best_params <- map(out, function(x) x[[idx]])

taus <- seq(-20, 40, 0.5)
p <- map_dbl(taus, function(t) nf_pdf(t, best_params$a, best_params$b, best_params$c, best_params$tmax))
sampled <- sample(taus, 1e4, replace = TRUE, prob = p)


inc_times <- rgamma(1e4, shape = params_real$inc_par2[["shape"]], scale = params_real$inc_par2[["scale"]])

si_posterior <- sampled + inc_times

invalid_si <- runif(1e4, -20, 40)

toss <- runif(1e4)

valid <- which(toss > best_params$pinvalid)

mixed <- c(
  si_posterior[valid], invalid_si[!valid]
)

p <- ggplot() +
  geom_density(aes(mixed, fill = "red"), alpha = 0.3, col = NA) +
  geom_histogram(aes(data_offset$si, y = ..density.., fill = "blue"), alpha = 0.3, col = NA) +
  scale_fill_identity(
    breaks = c("red", "blue"), labels = c("Posterior SI", "Data"),
    guide = "legend"
  ) +
  xlab("") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())


cowplot::save_plot("figures/cowling_nf_distr.png", p)