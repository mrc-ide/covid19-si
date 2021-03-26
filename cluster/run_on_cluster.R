library(context)
rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "context"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue")
source_files <- c("global2.R", "utils.R", "simulate2.R")
ctx <-context_save(
  root, packages = packages, sources = source_files
)
context::context_log_start()
misspec_offset <- -7
obj <- didehpc::queue_didehpc(ctx)
grp <- obj$lapply(
             index, function(i) {
               params_inc <- params_inc_all[[i]]
               params_offset <- params_offsets_all[[i]]
               params_iso <- params_iso_all[[i]]
               sim_data <- sampled[[i]]
               si_vec <- seq(params_offset + 0.5, max_valid_si, 1)
               width <- 0.1
               sim_data <- dplyr::arrange(sim_data, nu)
               fit_4a <- stan(
                 file = "scenario4a_mixture.stan",
                 data = list(
                   N = length(sim_data$si),
                   si = sim_data$si,
                   nu = sim_data$nu,
                   max_shed = max_shed,
                   offset1 = misspec_offset,
                   alpha2 = params_inc[["shape"]],
                   beta2 = 1 / params_inc[["scale"]],
                   alpha_invalid = alpha_invalid,
                   beta_invalid = beta_invalid,
                   max_valid_si = max_valid_si,
                   min_valid_si = params_offset,
                   min_invalid_si = min_invalid_si,
                   max_invalid_si = max_valid_si,
                   width = width,
                   M = length(si_vec),
                   si_vec = si_vec,
                   first_valid_nu = 1
                 ),
                 chains = 2, iter = 2000,
                 seed = 42,
                 verbose = TRUE
               )
               # outfile <- glue::glue("{prefix}_{index}.rds")
               # saveRDS(fit_4a, outfile)
               fit_4a
             }
           )
