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
# [ init:id   ]  c8d99ff52b3a84a33ed60bad1db84f35
# [ init:db   ]  rds
# [ init:path ]  relaxed_priors
# [ save:id   ]  2353ad978b96e51d71e527275d6cb0a0
# [ save:name ]  paleoclimatological_beardeddragon

obj <- didehpc::queue_didehpc(ctx)
fits <- obj$enqueue_bulk(model_features, fit_model)
tb1 <- obj$task_bundle_get('feminist_corydorascatfish')
fits <- tb1$results()
outfiles <- glue('stanfits/relaxed_priors/{model_features$model_prefix}_nf_fit.rds')
purrr::walk2(fits, outfiles, function(x, y) saveRDS(x, y))

# 01.07.2021 S3s4 mixture, only need run s4 models on this 'genial_stickinsect'
fits_s3s4 <- obj$enqueue_bulk(model_features[model_features$right_bias, ], s3s4_model)
# Fits with discrete pairs
root <- "pairs"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue", 'BH', 'RcppEigen')
source_files <- c("setup.R", 'cowling-data-prep.R')

ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = conan::conan_sources('./BH_1.75.0-0.tar.gz')
)
# [ init:id   ]  2244c4acbf66d27e2d3ca5dc8f1bdffa
# [ init:db   ]  rds
# [ init:path ]  pairs
# [ save:id   ]  72c43c2cd369394b24423a1cee2e2fe2
# [ save:name ]  draconic_finch
objpairs <- didehpc::queue_didehpc(ctx)
fits_pairs <- objpairs$enqueue_bulk(model_features, fit_model_to_pairs)
# 01.07.2021 Fit all models to discrete pairs only 'premythical_americankestrel'