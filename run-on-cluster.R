library(context)
rstan::rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "maxshed21_nf_priors"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue", 'BH', 'RcppEigen')
source_files <- c("setup.R", 'cowling-data-prep.R')

ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = conan::conan_sources('./BH_1.75.0-0.tar.gz')
)
# [ init:id   ]  3e274185e77c1fdff7be3052b1863a6d
# [ init:db   ]  rds
# [ init:path ]  maxshed21_nf_priors
# [ save:id   ]  ce5a01d1d2c918af36130fc3d93e4fe7
# [ save:name ]  unspirited_angwantibo

obj <- didehpc::queue_didehpc(ctx)
# 03.07.2021 Sensitivity analyses: max-shed 28 'nonscholastical_hound'
# 05.07.2021 max shed 21 and Neil's priors 'felicific_rhea'
fits <- obj$enqueue_bulk(model_features, fit_model)
tb1 <- obj$task_bundle_get('felicific_rhea')
fits <- tb1$results()
outfiles <- glue('stanfits/maxshed21_nfpriors/release/{model_features$model_prefix}_nf_fit.rds')
purrr::walk2(fits, outfiles, function(x, y) saveRDS(x, y))

# 01.07.2021 S3s4 mixture, only need run s4 models on this 'genial_stickinsect'
# 03.07.2021 Sensitivity analyses: max-shed 28 'bleareyed_garpike'
# 05.07.2021 max shed 21 and Neil's priors 'trashy_nilgai'
fits_s3s4 <- obj$enqueue_bulk(model_features[model_features$right_bias, ], s3s4_model)
fits_s3s4 <- obj$task_bundle_get('trashy_nilgai')
outfiles <- glue('stanfits/maxshed21_nfpriors/s3s4mix/{model_features$model_prefix[model_features$right_bias]}_nf_fit.rds')
purrr::walk2(fits_s3s4$results(), outfiles, function(x, y) saveRDS(x, y))

# Fits with discrete pairs
fits_pairs <- obj$enqueue_bulk(model_features, fit_model_to_pairs)
# 01.07.2021 Fit all models to discrete pairs only 'premythical_americankestrel'
# 03.07.2021 Sensitivity analyses: max-shed 28 'sporadic_goshawk'
# 05.07.2021 max shed 21 and Neil's priors 'zoometrical_xantusmurrelet'
# 08.07.2021 max shed 21 and Neil's priors s3s4mix, discrete pairs operose_anemonecrab
fits <- obj$task_bundle_get('operose_anemonecrab')
outfiles <- glue('stanfits/maxshed21_nfpriors/s3s4_discrete_pairs/{model_features$model_prefix}_nf_fit.rds')
purrr::walk2(fits$results(), outfiles, function(x, y) saveRDS(x, y))
