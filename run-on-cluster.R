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
# [ init:id   ]  11c4dfc686e0f0f77e65012ddebdb350
# [ init:db   ]  rds
# [ init:path ]  skew_normal
# [ save:id   ]  19bf1fa38fb7d6ff761c5d09ea1e7f0b
# [ save:name ]  beribboned_asiaticlesserfreshwaterclam
obj <- didehpc::queue_didehpc(ctx)
# 10.07.2021 'contemplative_arcticfox'
fits <- obj$enqueue_bulk(model_features, fit_model)
# 10.07.2021 'itchy_narwhale'
fits_pairs <- obj$enqueue_bulk(model_features, fit_model_to_pairs)
# 10.07.2021 s3s4 model 'combative_milksnake'
fits_s3s4 <- obj$enqueue_bulk(model_features[model_features$right_bias, ], s3s4_model)
# 10.07.2021 s3s4 model to discrete pairs 'cadaveric_cuscus'
fits_s3s4pairs <- obj$enqueue_bulk(model_features[model_features$right_bias, ], fit_model_to_s3s4pairs)
