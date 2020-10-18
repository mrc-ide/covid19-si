x <- -2
ulim <- x

inf_density <- dbeta(0.1/(max_shed + offset), shape1 = 2, shape2 = 6, log = TRUE) + log(0.1) - log(max_shed + offset)
inc_density <- dgamma((x - (-offset + 0.1)), shape = params_inc_og$shape, scale = params_inc_og$scale) + log(0.1)
out <-  exp(inf_density + inc_density)
s <- -offset + 0.2
while(s < ulim) {
print(s)
inf_density <- dbeta((s + offset)/(max_shed + offset), shape1 = 2, shape2 = 6, log = TRUE) + log(0.1) - log(max_shed + offset)
inc_density <- dgamma(x - s, shape = params_inc_og$shape, scale = params_inc_og$scale) + log(0.1)
print(inc_density)
out <- out + exp(inf_density + inc_density)
s <- s + 0.1
}
out <- log(out)
