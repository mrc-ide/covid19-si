fits <- pmap(model_features, fit_model)

## On the cluster
## fits <- obj$enqueue_bulk(model_features, fit_model)
