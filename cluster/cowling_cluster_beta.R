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

params_offset <- -11
params_inc <- params_real$inc_par2
max_shed <- params_real$maxshed3

####################
# read in the data #
####################

data <- readRDS("data/cowling_data_clean.rds")

real_data <- data %>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)%>%
  dplyr::filter(!is.na(nu))

real_data <- real_data[order(real_data$nu),]

saveRDS(real_data, "data/real_data_filtered.rds")

first_valid_nu <- match(params_offset + 2, real_data$nu) #changed to 2 for params_offset = 11 bc no nu = 10


######################
# fit the stan model #
######################
si_vec <- seq(params_offset + 1, max_si)
# real_data <- readRDS("data/real_data_filtered.rds")
# "13a461dbae54356ca011b022da92a96a"
beta_4a_recall <- obj$enqueue(stan(
  file = "scenario4arecall_mixture_beta.stan",
  data = list(
    N = length(real_data$si),
    si = real_data$si,
    nu = real_data$nu,
    max_shed = max_shed,
    offset1 = params_offset,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    max_invalid_si = max_si,
    min_invalid_si = min_invalid_si,
    width = width,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu =  first_valid_nu
  ),
  seed = 42,
  verbose = TRUE
))

beta_4a_recall$log()
# d6a42e2ec4d8fb9b31cf897902e30900
fit_4a_real_beta_11 <- obj$enqueue(stan(
  file = "scenario4a_mixture_beta.stan",
  data = list(
    N = length(real_data$si),
    si = real_data$si,
    nu = real_data$nu,
    max_shed = max_shed,
    offset1 = params_offset,
    alpha2 = params_inc[["shape"]],
    beta2 = 1 / params_inc[["scale"]],
    max_invalid_si = max_si,
    min_invalid_si = min_invalid_si,
    width = width,
    M = length(si_vec),
    si_vec = si_vec,
    first_valid_nu =  first_valid_nu
  ),
  seed = 42,
  verbose = TRUE
))
