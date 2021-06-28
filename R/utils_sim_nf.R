# simulate SI data using nf dist

# fitted parameters --> sampled TOST and SI distributions
sim_nf <- function(a, b, c, tmax, taus = seq(-20, 21, 0.1), n = 1e5,
                            shape_inc = params_real$inc_par2[["shape"]],
                            scale_inc = params_real$inc_par2[["scale"]],
                   shape_isol = 3,
                   scale_isol = 2,
                   offset_isol = 5){
  # TOST
  TOST <- rnf(n = n,
              taus = taus,
              a = a,
              b = b,
              c = c,
              tmax = tmax)
  
  # incubation period
  inc <- rgamma(shape = shape_inc,
                scale = scale_inc,
                n = n)
  
  # serial interval
  SI <- inc + TOST
  
  # delay onset to isolation (nu)
  nu <- rgamma(n = n, shape = shape_isol, scale = scale_isol) - offset_isol
  
  SI_full <- data.frame(SI = SI, TOST = TOST, inc = inc, nu = nu) 
  
  # filter on isolation
  SI_observed <- SI_full %>%
    filter(TOST < nu)
  
  # tags
  
  SI_full$tag <- "full"
  SI_observed$tag <- "obs"

  
  return(list(TOST = TOST,
              SI_full = SI_full,
              SI_observed = SI_observed))
  
}

## simulate the data (observed and underlying)

sim_si <- sim_nf(a = 5, b = 1, c = 0.2, tmax = 0, taus = seq(-20, 21, 0.1), n = 1e4,
                   shape_inc = params_real$inc_par2[["shape"]],
                   scale_inc = params_real$inc_par2[["scale"]],
                   shape_isol = 3,
                   scale_isol = 2,
                   offset_isol = 5)

SI_comp <- rbind(sim_si$SI_full, sim_si$SI_observed)

ggplot(SI_comp)+
  geom_histogram(aes(x = SI, y = ..density.., fill = tag), position = "identity" , binwidth = 1, alpha = 0.2)+
  theme_minimal()

saveRDS(file = "sim_data_manuscript.rds", sim_si)

