alpha1 <- 9.52381
beta1 <- 15.47619
alpha2 <- 9
beta2 <- 1 / 0.6666667
max_si <- 60

pdf_given_nu <- function(x, nu, max_shed, offset1, alpha1, beta1, alpha2,
                      beta2, width) {
  if (x > max_shed) ulim <- max_shed
  else ulim <- x
  if (ulim > nu) ulim <- nu
  max_shed_shifted <- max_shed - offset1;
  nu_shifted <- nu - offset1;
  f <- function(s) {
    dbeta((s - offset1) / max_shed_shifted, alpha1, beta1) *
      dgamma( (x - s), shape = alpha2, rate = beta2)
  }
  subd <- length(seq(offset1, ulim, width))
  out <- integrate(f, lower = offset1, upper = ulim, subdivisions = subd)
  out <- log(out$value)
  if(nu < max_shed) {
    out = out -
      pbeta(nu_shifted / max_shed_shifted, alpha1, beta1, log.p = TRUE)
  }
  out
}

normalisation_constant <- function(min_si, max_si, nu, max_shed,
                                   offset1, alpha1, beta1, alpha2,
                                   beta2, width) {
  si_vec <- seq(min_si, max_si, 1) + 0.5
  si_vec <- si_vec[si_vec <= max_si]
  out <- map_dbl(
    si_vec, function(x) pdf_given_nu(x, nu, max_shed,
                                     offset1, alpha1, beta1, alpha2,
                                     beta2, width)
  )
  exp(out)
}


normalisation_constant_stan <- function(min_si, max_si, nu, max_shed,
                                   offset1, alpha1, beta1, alpha2,
                                   beta2, width) {
  si_vec <- seq(min_si, max_si, 1) + 0.5
  si_vec <- si_vec[si_vec <= max_si]
  out <- map_dbl(
    si_vec, function(x) scenario4a_lpdf(x, nu, max_shed,
                                        offset1, alpha1, beta1, alpha2,
                                        beta2, width)
  )
  exp(out)
}




rstan::expose_stan_functions("stan-models/likelihoods.stan")
si_nu <- expand.grid(si = seq(-1, 60, 1) + 0.5, nu = 1:23)

si_nu$pdf_r <- pmap_dbl(
  si_nu, function(si, nu) {
    basic_pdf(si, nu, 21, -1, alpha1, beta1, alpha2, beta2, 0.1)
  }
)

si_nu$pdf_stan <- pmap_dbl(
  si_nu[, c("si", "nu")], function(si, nu) {
    scenario4a_lpdf(si, nu, 21, -1, alpha1, beta1, alpha2, beta2, 0.1)
  }
)

## PDF as calculated here and as calculated in Stan agree.
## goes almost to 0 for SI > 25, for all values of nu.
ggplot(si_nu) +
  geom_line(aes(si, exp(pdf_r), group = nu)) +
  geom_line(aes(si, exp(pdf_stan), group = nu), col = "red")




si_nu$den_stan <- map_dbl(
  si_nu$nu, function(nu) {
    s4_normalising_constant(nu, 21, -1, alpha1, beta1, alpha2, beta2,
                            max_si, 0.1)
  }
)

si_nu$den_r <- map_dbl(
  si_nu$nu, function(nu) {
    sum(normalisation_constant(-1, 20, nu, 21, -1, alpha1, beta1, alpha2,
                           beta2, 0.1))
  }
)
## Normalised density; almost the same for all values of nu
si_nu$normalised_r <- si_nu$pdf_r - log(si_nu$den_r)
si_nu$normalised_stan <- si_nu$pdf_stan - log(si_nu$den_stan)


ggplot(si_nu) +
  geom_line(aes(si, normalised_r, group = nu)) +
  geom_line(aes(si, normalised_stan, group = nu), col = "red")


si_nu[si_nu$nu %in% 11, ] %>%
ggplot() +
  geom_line(aes(si, normalised_r, group = nu)) +
  geom_line(aes(si, pdf_r, group = nu), linetype = "dashed")




outfiles <- "data/4a_mix_with_normalisation_sim_1data.rds"
mixed <- readRDS(outfiles)
mixed <- mutate_if(mixed, is.numeric, round)
mixed <- mixed[mixed$si != 0, ]
mixed <- mixed[mixed$nu != 0, ]
## For now let min_invalid_si be -1
##mixed <- mixed[mixed$si > -1, ]
## grid likelihood
grid <- expand.grid(
  a1 = 5:15, b1 = 5:20, pinvalid = seq(0, 1, by = 0.1)
)

grid$likelihood <- pmap_dbl(
  grid,
  function(a1, b1, pinvalid) {
    out <- pmap_dbl(
      mixed[, c("si", "nu")],
      function(si, nu) {
        valid <- pdf_given_nu(si, nu, 21, -1, a1, b1, alpha2,
                              beta2, 0.1)
        den <- sum(normalisation_constant(-1, max_si, nu, 21, -1,
                                          a1, b1, alpha2,
                                          beta2, 0.1))
        valid <- valid - log(den)
        invalid <- log(pinvalid) - log(max_si + 1)
        invalid + log(1 - pinvalid) + valid
      }
    )
    sum(out)
  }
)


grid <- expand.grid(
  a1 = 1:15, b1 = beta1, pinvalid = seq(0, 1, by = 0.01)
)

grid$likelihood <- pmap_dbl(
  grid,
  function(a1, b1, pinvalid) {
    out <- pmap_dbl(
      mixed[, c("si", "nu")],
      function(si, nu) {
        invalid <-  log(max_si + 1)

        if (si <= -1) {
          return(log(pinvalid) - invalid)
        } else {
          valid <- pdf_given_nu(si, nu, 21, -1, a1, b1, alpha2,
                                beta2, 0.1)
          den <- sum(normalisation_constant(-1, max_si, nu, 21, -1,
                                            a1, b1, alpha2,
                                            beta2, 0.1))
          valid <- valid - log(den)
        }
        log(pinvalid) - invalid + log(1 - pinvalid) + valid
      }
    )
    sum(out)
  }
)


grid_finite <- grid[is.finite(grid$likelihood), ]

ggplot(grid_finite, aes(pinvalid, likelihood)) +
  geom_line()

ggplot(grid_finite, aes(x = a1, y = b1, z = ))



############
### Check that we get the same answer as from the previous
### implementation, which we believe to be correct
si_vec <- seq(-0.5, 25, 0.5)
nu_vec <- 0:30

out <- pdf_matrix(nu_vec, si_vec, 21, -1, alpha1, beta1, alpha2, beta2, 0.1)
current <- colSums(out)

previous <- map_dbl(nu_vec, function(nu) {
  s4_normalising_constant(nu, 21, -1, alpha1, beta1, alpha2, beta2, 25, 0.1)
})


## Test with some values of nu repeated
nu_vec <- c(5, 5, 7, 7, 7, 8, 9)

out <- pdf_matrix(nu_vec, si_vec, 21, -1, alpha1, beta1, alpha2, beta2, 0.1)
current <- colSums(out)
