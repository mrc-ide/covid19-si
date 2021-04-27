

#1. read in rds files

fit3 <- readRDS("stanfits/release/scenario3a_mixture_nf.rds")
fit3_recall <- readRDS("stanfits/release/scenario3arecall_mixture_nf.rds")
fit4 <- readRDS("stanfits/release/scenario4a_mixture_nf.rds")
fit4_recall <- readRDS("stanfits/release/scenario4arecall_mixture_nf.rds")


## processing fits --> tables and figures

# testing

tab1_s3 <- table1_fun(fit = fit3)
tab1_s4 <- table1_fun(fit = fit4)

samples_s3 <- sample_dist_fun(tab1_s3,
                              fit = rstan::extract(fit3))
samples_s4 <- sample_dist_fun(tab1_s4,
                              fit = rstan::extract(fit4))
tab2_s3 <- table2_fun(samples_s3)
tab2_s4 <- table2_fun(samples_s4)

TOST3 <- samples_s3$TOST_bestpars
TOST4 <- samples_s4$TOST_bestpars

p <- TOST_fig_fun(TOST3)

SI3 <- samples_s3$SI_bestpars$SI
SI4 <- samples_s4$SI_bestpars$SI
data <- readRDS("data/cowling_data_clean.rds")
data <- data%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)%>%
  dplyr::filter(!is.na(nu))

SI4_con <- pred_observed_SI_fun(SI = samples_s4$SI_bestpars,
                                data = data,
                                mixture = TRUE,
                                recall = FALSE,
                                isol = TRUE,
                                tab1 = tab1_s4,
                                n = 1e5)
SIcomp_fig_fun(SI1 = SI4, SI2 = SI4_con, data = data)
TOST3_post <- samples_s3$TOST_post
TOST3_post <- reshape2::melt(TOST3_post[,1:1000])
p <- ggplot(TOST3_post)
  

 p <- ggplot()+
  geom_density(data = TOST3_post, aes(x = value, y = ..density.., group = variable), size = 1, colour = "grey", alpha = 0.3)+
  geom_density(aes(x = TOST3, y = ..density..), colour = "black")+
  theme_minimal()+
  xlab("TOST (days)")
 
 q <- ggplot()+
   stat_ecdf(data = TOST3_post, aes(x = value, group = variable), size = 1, colour = "grey", alpha = 0.3)+
   stat_ecdf(aes(x = TOST3), colour = "black")+
   theme_minimal()+
   xlab("TOST (days)")

