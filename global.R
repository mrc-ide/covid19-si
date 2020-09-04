library(dplyr)
library(rstan)
library(drake)
library(hermione)
library(ggmcmc)
library(ggplot2)
library(epitrix)
library(purrr)
max_shed <- 21

## https://doi.org/10.1016/S1473-3099(20)30287-5
mean_inf <- 6.3 ## days
sd_inf <- 4.2 ## days
params_inf <- epitrix::gamma_mucv2shapescale(mean_inf, sd_inf/ mean_inf)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7014672/
mean_inc <- 6.5 ## days
sd_inc <- 2.6 ## days
params_inc <- epitrix::gamma_mucv2shapescale(mean_inc, sd_inc/ mean_inc)