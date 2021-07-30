library(context)
rstan::rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "gamma"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix",
              "glue", 'BH', 'RcppEigen')
source_files <- c("setup.R", 'cowling-data-prep.R')

ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = conan::conan_sources('./BH_1.75.0-0.tar.gz')
)
# [ init:id   ]  bf6f92ea59df96c383cfe0a629dbbd3a
# [ init:db   ]  rds
# [ init:path ]  gamma
# [ save:id   ]  d997cee68365a60ea8afc0eaf48bc198
# [ save:name ]  polluted_mangabey
obj <- didehpc::queue_didehpc(ctx)
# 27.07.2021 'sudden_barasinga'
# fits <- obj$enqueue_bulk(model_features, fit_model)
fits <- obj$task_bundle_get('sudden_barasinga')
# recall_3 'purgatorial_lovebird'
recall_3 <-  obj$enqueue_bulk(model_features[6, ], fit_model)
# 27.07.2021 'tridymite_geese'
# fits_pairs <- obj$enqueue_bulk(model_features, fit_model_to_pairs)
fits_pairs <- obj$task_bundle_get('tridymite_geese')
# recall3_pairs 'nasty_gull'
recall3_pairs <- obj$enqueue_bulk(model_features[6, ], fit_model_to_pairs)
# 27.07.2021 s3s4 model 'gullible_rasbora'
# fits_s3s4 <- obj$enqueue_bulk(model_features[model_features$right_bias, ], s3s4_model)
fits_s3s4 <- obj$task_bundle_get('gullible_rasbora')
# 27.07.2021 s3s4 model to discrete pairs 'pseudoliterary_tomtit'
# fits_s3s4pairs <- obj$enqueue_bulk(model_features[model_features$right_bias, ], fit_model_to_s3s4pairs)
fits_s3s4pairs <- obj$task_bundle_get('pseudoliterary_tomtit')

outfiles <- glue('stanfits/gamma/s3s4pairs/{model_features$model_prefix[model_features$right_bias]}_gamma_fit.rds')
purrr::walk2(fits_s3s4pairs$results(), outfiles, function(x, y) saveRDS(x, y))


outfiles <- glue('stanfits/gamma/discrete_pairs/{model_features$model_prefix}_gamma_fit.rds')
purrr::walk2(fits_pairs$results(), outfiles, function(x, y) saveRDS(x, y))

pwalk(
  model_features,
  function() {
    outfile <- glue('stanfits/skew_normal/discrete_pairs/{model_features$model_prefix}_skew_normal_fit.rds')
    out <- fit_model_to_pairs
    saveRDS(out, outfile)
  }
)
