library(context)
rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "4a_mixture"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue")
source_files <- c("global2.R", "utils.R", "scenario4a_cluster.R")
ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = provisionr::package_sources(local = "BH_1.75.0-0.zip")
  )

context::context_log_start()

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
                   offset1 = params_offset,
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
               fit_4a
             }
           )

tb <- obj$task_bundle_get('pseudomedical_swallowtailbutterfly')
t <- obj$task_get('251f7692e0fa235fa827a8c46f2d6c70')
