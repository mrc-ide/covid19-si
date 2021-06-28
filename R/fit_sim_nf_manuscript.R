## fitting s3 and s4 to simulated data with isolation

## read in simulated data
sim_si <- readRDS("sim_data_manuscript.rds")
si_vec <- seq(-20, 40, 1)

## Make sure data are arranged in order of nu
data <- sim_si$SI_observed
data <- arrange(data, nu)

#fit s3
fit3_sim <- stan(
  file = here::here("stan-models/scenario3a_mixture_nf.stan"),
  data = list(
    N = nrow(data),
    si = data$SI,
    nu = data$nu,
    max_shed = 21,
    alpha2 = params_real$inc_par2[["shape"]],
    beta2 = 1 / params_real$inc_par2[["scale"]],
    max_invalid_si = 40,
    min_invalid_si = -20,
    width = 1,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
    ##tmax = 0
  ),
  chains = 2, iter = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

tab1_s3 <- table1_fun(fit3_sim)
best_s3 <- tab1_s3$best
names(best_s3) <- rownames(tab1_s3)

fit_s3_sim <- sim_nf(a = best_s3["a"], b = best_s3["b"], c = best_s3["c"], tmax = best_s3["tmax"], taus = seq(-20, 21, 0.1), n = 1e5,
                 shape_inc = params_real$inc_par2[["shape"]],
                 scale_inc = params_real$inc_par2[["scale"]],
                 shape_isol = 3,
                 scale_isol = 2,
                 offset_isol = 5) # don't worry about iso - just plot full SI

ggplot(SI_comp)+
  geom_histogram(aes(x = SI, y = ..density.., fill = tag), position = "identity" , binwidth = 1, alpha = 0.2)+
  geom_density(data = fit_s3_sim$SI_full, aes(x = SI), size = 1, colour = "grey")+
  theme_minimal()


## fit s4
fits_4a <- stan(
  file = here::here("stan-models/scenario4arecall_mixture_nf.stan"),
  data = list(
    N = nrow(data),
    si = data$SI,
    nu = data$nu,
    max_shed = 21,
    alpha2 = params_real$inc_par2[["shape"]],
    beta2 = 1 / params_real$inc_par2[["scale"]],
    max_invalid_si = 40,
    min_invalid_si = -20,
    width = 1,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
    ##tmax = 0
  ),
  chains = 2, iter = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


## plot comparisons

