

## read in rds files

fit3 <- readRDS("stanfits/release/scenario3a_mixture_nf.rds")
fit3_recall <- readRDS("stanfits/release/scenario3arecall_mixture_nf.rds")
fit4 <- readRDS("stanfits/release/scenario4a_mixture_nf.rds")
fit4_recall <- readRDS("stanfits/release/scenario4arecall_mixture_nf.rds")
fit3_nomix_recall <- readRDS("stanfits/release/scenario3arecall_no_mixture_nf.rds")
fit4_nomix_recall <- readRDS("stanfits/release/scenario4arecall_no_mixture_nf.rds")
fit3_leftbias <- readRDS("stanfits/release/scenario3a_mixture_left_bias_nf.rds")
fit3_nomix_leftbias <- readRDS("stanfits/release/scenario3a_nomixture_left_bias_nf.stan.rds")

## processing fits --> tables and figures

 tab1_s3mix <- table1_fun(fit3)
 tab1_s4mix <- table1_fun(fit4)
 tab1_s3mixrecall <- table1_fun(fit3_recall)
 tab1_s4mixrecall <- table1_fun(fit4_recall)
 tab1_s3recall <- table1_fun(fit3_nomix_recall)
 tab1_s4recall <- table1_fun(fit4_nomix_recall)
 tab1_s3mixleftbias <- table1_fun(fit3_leftbias)
 tab1_s3leftbias <- table1_fun(fit3_nomix_leftbias)

## sample distributions

 samples_s3mix <- sample_dist_fun(tab1_s3mix, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit3),
                             shape_inc = params_real$inc_par2[["shape"]],
                             scale_inc = params_real$inc_par2[["scale"]])
 samples_s4mix <- sample_dist_fun(tab1_s4mix, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit4),
                                  shape_inc = params_real$inc_par2[["shape"]],
                                  scale_inc = params_real$inc_par2[["scale"]])
 samples_s3mixrecall <- sample_dist_fun(tab1_s3mixrecall, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit3_recall),
                                  shape_inc = params_real$inc_par2[["shape"]],
                                  scale_inc = params_real$inc_par2[["scale"]])
 samples_s4mixrecall <- sample_dist_fun(tab1_s4mixrecall, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit4_recall),
                                  shape_inc = params_real$inc_par2[["shape"]],
                                  scale_inc = params_real$inc_par2[["scale"]])
 samples_s3recall <- sample_dist_fun(tab1_s3recall, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit3_nomix_recall),
                                  shape_inc = params_real$inc_par2[["shape"]],
                                  scale_inc = params_real$inc_par2[["scale"]])
 samples_s4recall <- sample_dist_fun(tab1_s4recall, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit4_nomix_recall),
                                  shape_inc = params_real$inc_par2[["shape"]],
                                  scale_inc = params_real$inc_par2[["scale"]])
 samples_s3mixleftbias <- sample_dist_fun(tab1_s3mixleftbias, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit3_leftbias),
                                        shape_inc = params_real$inc_par2[["shape"]],
                                        scale_inc = params_real$inc_par2[["scale"]])
 samples_s3leftbias <- sample_dist_fun(tab1_s3leftbias, taus = seq(-20, 40, 0.1), n = 1e4, rstan::extract(fit3_nomix_leftbias),
                                          shape_inc = params_real$inc_par2[["shape"]],
                                          scale_inc = params_real$inc_par2[["scale"]])

## table 2 - summary stats for sampled distributions

 tab2_s3mix <- table2_fun(samples_s3mix)
 tab2_s4mix <- table2_fun(samples_s4mix)
 tab2_s3mixrecall <- table2_fun(samples_s3mixrecall)
 tab2_s4mixrecall <- table2_fun(samples_s4mixrecall)
 tab2_s3recall <- table2_fun(samples_s3recall)
 tab2_s4recall <- table2_fun(samples_s4recall)
 tab2_s3mixleftbias <- table2_fun(samples_s3mixleftbias)
 tab2_s3leftbias <- table2_fun(samples_s3leftbias)


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

 SI4_con <- expected_SI_fun(SI = samples_s4$SI_bestpars,
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


