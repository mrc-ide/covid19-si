library(context)
rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "4a_mixture"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue")
source_files <- "global2.R"
ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = provisionr::package_sources(local = "BH_1.75.0-0.zip")
)

context::context_log_start()

obj <- didehpc::queue_didehpc(ctx)


data_offset <- readRDS('cowling_data_cleaned.rds')
# fit the model
si_vec <- seq(-20, 40, 1)
## Scenario3 NF Distribution
t3 <- obj$enqueue(stan(
  file = "scenario3a_mixture_nf.stan",
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    nu = data_offset$nu,
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
  seed = 42,
  verbose = TRUE
))





## Scenario3 + recall NF Distribution
fit_3arecall <- obj$enqueue(stan(
  file = "scenario3arecall_mixture_nf.stan",
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    nu = data_offset$nu,
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
  ), verbose = TRUE
))



### Scenario 4  NF Distribution
fit_4a <- obj$enqueue(stan(
  file = "scenario4a_mixture_nf.stan",
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    nu = data_offset$nu,
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
  ), verbose = TRUE
))

### Scenario 4 + recall  NF Distribution
fit_4arecall <- obj$enqueue(stan(
  file = "scenario4arecall_mixture_nf.stan",
  data = list(
    N = nrow(data_offset),
    si = data_offset$si,
    nu = data_offset$nu,
    max_shed = 21,
    alpha2 = params_real$inc_par2[["shape"]],
    beta2 = 1 / params_real$inc_par2[["scale"]],
    max_invalid_si = 40,
    min_invalid_si = -20,
    width = 1,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu = 1
  ), verbose = TRUE
))
