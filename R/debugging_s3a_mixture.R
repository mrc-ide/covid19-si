# debugging S3a + mixture

# first lets look at the likelihoods with and without widths
expose_stan_functions("stan-models/likelihoods.stan")

SI <- c(-11:24)
N <- length(SI)
out_nowidth <- vector(length = N)
out_width <- vector(length = N)

# without width window multiplier
for(i in 1:N){
 out_nowidth[i] <- scenario3a_nowidth_lpdf(x = SI[i], 
                        max_shed = 21, 
                        offset = 3,
                        alpha1 = 2, 
                        beta1 = 6, 
                        alpha2 = params_inc_og$shape, 
                        beta2 = params_inc_og$scale, 
                        width = 0.1)

}

# with width window multiplier
for(i in 1:N){
  out_width[i] <- scenario3a_width_lpdf(x = SI[i], 
                                            max_shed = 21, 
                                            offset = 3,
                                            alpha1 = 2, 
                                            beta1 = 6, 
                                            alpha2 = params_inc_og$shape, 
                                            beta2 = params_inc_og$scale, 
                                            width = 0.1)
  
}
# expected difference
width <- 0.1
max_shed <- 21
offset <- 3

expected_diff <- log(width/(max_shed + offset)) + log(width)

# difference

out_diff <- out_width - out_nowidth


# now lets look at log likelihood of being invalid for the same SIs
out_invalid <- vector(length = N)

for(i in 1:N){
out_invalid[i] <- invalid_lpdf(x = SI[i], max_si = max(SI) + 0.001, min_si = min(SI) - 0.001, alpha_invalid = 0.5,
             beta_invalid = 0.5)
}

# bringing this together in the mixture model
target_width <- vector(length = N)
p_inv <- 0.99
for (i in 1:N) {
  if(SI[i]> -offset){
    target_width[i] <- p_inv*exp(out_invalid[i]) + (1-p_inv)*exp(out_width[i])
  } else{
    target_width[i] <- p_inv*exp(out_invalid[i])
  }
}

plot(target_width)
plot(log(target_width))

target_nowidth <- vector(length = N)

for (i in 1:N) {
  if(SI[i]> -offset){
    target_nowidth[i] <- p_inv*exp(out_invalid[i]) + (1-p_inv)*exp(out_nowidth[i])
  } else{
    target_nowidth[i] <- p_inv*exp(out_invalid[i])
  }
}
  
plot(target_nowidth)
plot(log(target_nowidth))

mat <- data.frame(log(target_width), log(target_nowidth))

matplot(mat, type = "l", xlab = "SI (days)", ylab = "loglik")

target_diff <- log(target_width) - log(target_nowidth) + out_invalid
