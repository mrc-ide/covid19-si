library(context)
rstan_options(auto_write = FALSE)
# config <- didehpc::didehpc_config(cluster = 'fi--didemrchnb')
options(didehpc.cluster = 'fi--didemrchnb')
root <- "simulated"
packages <- c("rstan", "dplyr","purrr", "ggplot2", "epitrix", "glue", 'RcppEigen')
source_files <- c("R/utils_process_fits.R", "R/utils.R", "R/utils_sim_nf.R")
src <- conan::conan_sources('local::BH_1.75.0-0.tar.gz')
ctx <-context_save(
  root, packages = packages, sources = source_files,
  package_sources = src
)
# [ open:db   ]  rds
# [ save:id   ]  32d2d01ee69329404984186a17bcdf2f
# [ save:name ]  atrophied_johndory
obj <- didehpc::queue_didehpc(ctx)




## tab1_s3 <- table1_fun(fit3_sim)
## best_s3 <- tab1_s3$best
## names(best_s3) <- rownames(tab1_s3)

## fit_s3_sim <- sim_nf(
##   a = best_s3["a"], b = best_s3["b"], c = best_s3["c"], tmax = best_s3["tmax"], taus = seq(-20, 21, 0.1), n = 1e5,
##   shape_inc = params_real$inc_par2[["shape"]],
##   scale_inc = params_real$inc_par2[["scale"]],
##   shape_isol = 3, scale_isol = 2, offset_isol = 5
## ) # don't worry about iso - just plot full SI

## ggplot(SI_comp)+
##   geom_histogram(aes(x = SI, y = ..density.., fill = tag), position = "identity" , binwidth = 1, alpha = 0.2)+
##   geom_density(data = fit_s3_sim$SI_full, aes(x = SI), size = 1, colour = "grey")+
##   theme_minimal()



task <- obj$enqueue(fit3_sim()) 
# 200f04111e6aabba82ba1b029b638bea


t4 <- obj$enqueue(fits_4a()) 
# 05e00ed215ec184bbe7eafe29aeca325
