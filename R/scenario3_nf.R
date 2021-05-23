fits <- pmap(model_features, fit_model)

## cowling data is now s3/s4 mix data. This will now run s3/s4
fits <- pmap(model_features[model_features$right_bias, ], fit_model)
## s3s4fits <- pmap(model_features[model_features$right_bias, ], fit_model)
## On the cluster
## fits <- obj$enqueue_bulk(model_features, fit_model)

##
leaky_fits <- pmap(model_features[model_features$right_bias, ], fit_leaky_model)
## leaky_fits <- obj$enqueue_bulk(model_features[model_features$right_bias, ], fit_leaky_model)





