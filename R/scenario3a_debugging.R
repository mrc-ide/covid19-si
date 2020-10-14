grid <- expand.grid(
  shape1 = seq(1,30, by = 0.5),
  shape2 = seq(1,30, by = 0.5)
)

funcs <- expose_stan_functions("stan-models/scenario3a.stan")
out <- vector(length = grid$shape1)
x <- seq(-2,24)
#x <- data_offset$si

for(j in 1:length(x)){

for(i in 1:length(grid$shape1)){
 out[i] <- scenario3a_lpdf(x = x[j], max_shed = 21, offset = 3, alpha2 = shape_inc, beta2 = rate_inc,
                alpha1 = grid$shape1[i], beta1 = grid$shape2[i]) 
}

summ <- sum(out == -Inf)
print(summ)
}
