library(context)
rstan::rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "skew_normal"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix",
              "glue", 'BH', 'RcppEigen')
source_files <- c("setup.R", 'cowling-data-prep.R')

ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = conan::conan_sources('./BH_1.75.0-0.tar.gz')
)
obj <- didehpc::queue_didehpc(ctx)
fits <- obj$enqueue_bulk(model_features, fit_model)
fits_pairs <- obj$enqueue_bulk(model_features, fit_model_to_pairs)
fits_s3s4 <- obj$enqueue_bulk(model_features[model_features$right_bias, ], s3s4_model)
fits_s3s4pairs <- obj$enqueue_bulk(model_features[model_features$right_bias, ], fit_model_to_s3s4pairs)
