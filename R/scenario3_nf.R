params_inc <- params_real$inc_par2
si_vec <- seq(-20, max_valid_si)
standata <- list(
  N = nrow(cowling_data), si = cowling_data$si, max_shed = max_shed,
  alpha2 = params_inc[["shape"]], beta2 = 1 / params_inc[["scale"]],
  M = length(si_vec), si_vec = si_vec, width = 0.1
)

fits <- pmap(
  model_features, function(mixture, recall, right_bias, model_prefix) {
    prefix <- glue("{model_prefix}_nf")
    infile <- glue("stan-models/{prefix}.stan")
    if (!file.exists(infile)) message("Does not exist ", infile)
    if(mixture) {
      standata$max_invalid_si <- max_invalid_si
      standata$min_invalid_si <- min_invalid_si
    }
    if (right_bias) standata$nu <- cowling_data$nu

    fit <- stan(
      file = infile, data = standata,  verbose = TRUE, iter = iter
    )

    fit
  }
)
