library(context)
rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "simulated"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue")
source_files <- c("R/utils_process_fits.R", "utils.R", "R/utils_sim_nf.R")
ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = provisionr::package_sources(local = "BH_1.75.0-0.zip")
)
context::context_log_start()
misspec_offset <- -7
obj <- didehpc::queue_didehpc(ctx)

#fit s3
fit3_sim <- function() {
 stan(
  file = here::here("stan-models/scenario3a_mixture_nf.stan"),
  data = list(
    N = nrow(sim_si),
    si = sim_si$SI,
    nu = sim_si$nu,
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
}

## tab1_s3 <- table1_fun(fit3_sim)
## best_s3 <- tab1_s3$best
## names(best_s3) <- rownames(tab1_s3)

## fit_s3_sim <- sim_nf(
##   a = best_s3["a"], b = best_s3["b"], c = best_s3["c"], tmax = best_s3["tmax"], taus = seq(-20, 21, 0.1), n = 1e5,
##   shape_inc = params_real$inc_par2[["shape"]],
##   scale_inc = params_real$inc_par2[["scale"]],
##   shape_isol = 3, scale_isol = 2, offset_isol = 5
## ) # don't worry about iso - just plot full SI

## ggplot(SI_comp)+
##   geom_histogram(aes(x = SI, y = ..density.., fill = tag), position = "identity" , binwidth = 1, alpha = 0.2)+
##   geom_density(data = fit_s3_sim$SI_full, aes(x = SI), size = 1, colour = "grey")+
##   theme_minimal()


## fit s4
fits_4a <- function() {
  stan(
  file = here::here("stan-models/scenario4arecall_mixture_nf.stan"),
  data = list(
    N = nrow(sim_si),
    si = sim_si$SI,
    nu = sim_si$nu,
    max_shed = 21,
    alpha2 = params_real$inc_par2[["shape"]],
    beta2 = 1 / params_real$inc_par2[["scale"]],
    max_invalid_si = 40,
    min_invalid_si = -20,
    width = 1,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
  ),
  chains = 2, iter = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
  )
}
