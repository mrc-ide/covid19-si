library(context)
rstan::rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "maxshed28_relaxed_priors"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue", 'BH', 'RcppEigen')
source_files <- c("setup.R", 'cowling-data-prep.R')

ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = conan::conan_sources('./BH_1.75.0-0.tar.gz')
)
# [ init:id   ]  36b8556ef3951ac6a93886240dcb9ce1
# [ init:db   ]  rds
# [ init:path ]  maxshed28_relaxed_priors
# [ save:id   ]  6eb2ad43578964f3d310bcab02373987
# [ save:name ]  unmystified_rainbowlorikeet

obj <- didehpc::queue_didehpc(ctx)
# 03.07.2021 Sensitivity analyses: max-shed 28 'nonscholastical_hound'
fits <- obj$enqueue_bulk(model_features, fit_model)
tb1 <- obj$task_bundle_get('nonscholastical_hound')
fits <- tb1$results()
outfiles <- glue('stanfits/relaxed_priors/{model_features$model_prefix}_nf_fit.rds')
purrr::walk2(fits, outfiles, function(x, y) saveRDS(x, y))

# 01.07.2021 S3s4 mixture, only need run s4 models on this 'genial_stickinsect'
# 03.07.2021 Sensitivity analyses: max-shed 28 'bleareyed_garpike'
fits_s3s4 <- obj$enqueue_bulk(model_features[model_features$right_bias, ], s3s4_model)
fits_s3s4 <- obj$task_bundle_get('bleareyed_garpike')
outfiles <- glue('stanfits/s3s4mix/{model_features$model_prefix[model_features$right_bias]}_nf_fit.rds')
purrr::walk2(fits_s3s4$results(), outfiles, function(x, y) saveRDS(x, y))

# Fits with discrete pairs
fits_pairs <- obj$enqueue_bulk(model_features, fit_model_to_pairs)
# 01.07.2021 Fit all models to discrete pairs only 'premythical_americankestrel'
# 03.07.2021 Sensitivity analyses: max-shed 28 'sporadic_goshawk'
fits <- objpairs$task_bundle_get('sporadic_goshawk')
outfiles <- glue('stanfits/discrete_pairs/{model_features$model_prefix}_nf_fit.rds')
purrr::walk2(fits$results(), outfiles, function(x, y) saveRDS(x, y))
