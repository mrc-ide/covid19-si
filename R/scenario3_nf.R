fits <- pmap(model_features, fit_model)

## cowling data is now s3/s4 mix data. This will now run s3/s4
fits <- pmap(model_features[model_features$right_bias, ], s3s4_model)
## s3s4fits <- pmap(model_features[model_features$right_bias, ], s3s4_model)
## On the cluster
## fits <- obj$enqueue_bulk(model_features, s3s4_model)

##
leaky_fits <- pmap(model_features[4, ], fit_leaky_model)
## leaky_fits <- obj$enqueue_bulk(model_features[model_features$right_bias, ], fit_leaky_model)


s3data$max_invalid_si <- max_invalid_si
s3data$min_invalid_si <- min_invalid_si

fit <- stan(
  file = "stan-models/scenario3a_mixture_skew_normal.stan",
  data = s3data, verbose = TRUE, chains = 2, iter = 5000
)

saveRDS(fit, "stanfits/skew_normal/scenario3a_mixture_skew_normal_fit.stan")

