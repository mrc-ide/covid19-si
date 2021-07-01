library(context)
rstan::rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "relaxed_priors"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue", 'BH', 'RcppEigen')
source_files <- c("setup.R", 'cowling-data-prep.R')

ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = conan::conan_sources('./BH_1.75.0-0.tar.gz')
)

obj <- didehpc::queue_didehpc(ctx)
fits <- obj$enqueue_bulk(model_features, fit_model)
tb1 <- obj$task_bundle_get('churnable_nightheron')
fits <- tb1$results()
outfiles <- glue('stanfits/relaxed_priors/{model_features$model_prefix}_nf_fit.rds')
purrr::walk2(fits, outfiles, function(x, y) saveRDS(x, y))
