

## read in rds files

fit3 <- readRDS("stanfits/release/scenario3a_mixture_nf_long_max_shed.rds")
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


## table 2 - summary stats for sampled distributions

 tab2_s3mix <- table2_fun(samples_s3mix)
 tab2_s4mix <- table2_fun(samples_s4mix)
 tab2_s3mixrecall <- table2_fun(samples_s3mixrecall)
 tab2_s4mixrecall <- table2_fun(samples_s4mixrecall)
 tab2_s3recall <- table2_fun(samples_s3recall)
 tab2_s4recall <- table2_fun(samples_s4recall)


## TOST plot
 TOST3mix <- samples_s3mix$TOST_bestpars
 TOST4mix <- samples_s4mix$TOST_bestpars
 TOST3recall <- samples_s3recall$TOST_bestpars
 TOST4recall <- samples_s4recall$TOST_bestpars
 TOST4mix_recall <- samples_s4mixrecall$TOST_bestpars
 TOST3mix_leftbias <- samples_s3mixleftbias$TOST_bestpars
 TOST3leftbias <- samples_s3leftbias$TOST_bestpars

 p1 <- TOST_fig_fun(TOST3mix)
 p2 <- TOST_fig_fun(TOST4mix)
 p3 <- TOST_fig_fun(TOST3recall)
 p4 <- TOST_fig_fun(TOST4recall)
 p5 <- TOST_fig_fun(TOST4mix_recall)
 p6 <- TOST_fig_fun(TOST3mix_leftbias)
 p7 <- TOST_fig_fun(TOST3leftbias)

## SI plot
 SI3mix <- samples_s3mix$SI_bestpars$SI
 SI4mix <- samples_s4mix$SI_bestpars$SI
 SI3recall <- samples_s3recall$SI_bestpars$SI
 SI4recall <- samples_s4recall$SI_bestpars$SI
 SI4mix_recall <- samples_s4mixrecall$SI_bestpars$SI
 SI3mix_leftbias <- samples_s3mixleftbias$SI_bestpars$SI
 SI3leftbias <- samples_s3leftbias$SI_bestpars$SI

 library(fitdistrplus)

 data <- readRDS("data/cowling_data_clean.rds")
 data <- data%>%
    mutate(si = as.numeric(si))%>%
    dplyr::rename(nu = onset_first_iso)%>%
    dplyr::filter(!is.na(nu))

 SI3mix_con <- expected_SI(SI = samples_s3mix$SI_bestpars,
                               data = data,
                               mixture = TRUE,
                               recall = FALSE,
                               isol = FALSE,
                               tab1 =  tab1_s3mix ,
                               n = 1e4
                               )
 SI4mix_con <- expected_SI(SI = samples_s4mix$SI_bestpars,
                                    data = data,
                                    mixture = TRUE,
                                    recall = FALSE,
                                    isol = TRUE,
                                    tab1 =  tab1_s4mix ,
                                    n = 1e4)
 SI4mix_con_nu0 <- expected_SI(SI = samples_s4mix$SI_bestpars,
                               data = data,
                               mixture = TRUE,
                               recall = FALSE,
                               isol = TRUE,
                               tab1 =  tab1_s4mix ,
                               n = 1e4)
 SI4mix_con_empiricnu <- expected_SI_empiricnu(SI = samples_s4mix$SI_bestpars,
                               data = data,
                               mixture = TRUE,
                               recall = FALSE,
                               isol = TRUE,
                               tab1 =  tab1_s4mix ,
                               n = 1e4)

 SI3recall_con <- expected_SI(SI = samples_s3recall$SI_bestpars,
                                       data = data,
                                       mixture = FALSE,
                                       recall = TRUE,
                                       isol = FALSE,
                                       tab1 =  tab1_s3recall,
                                       n = 1e4)
 SI4recall_con <- expected_SI(SI = samples_s4recall$SI_bestpars,
                                       data = data,
                                       mixture = FALSE,
                                       recall = TRUE,
                                       isol = TRUE,
                                       tab1 = tab1_s4recall,
                                       n = 1e4)
 SI4mix_recall_con <- expected_SI(SI = samples_s4mixrecall$SI_bestpars,
                                           data = data,
                                           mixture = TRUE,
                                           recall = TRUE,
                                           isol = TRUE,
                                           tab1 = tab1_s4mixrecall,
                                           n = 1e4)

 p1 <- SIcomp_fig_fun(SI1 = SI3mix, SI2 = SI3mix, data = data)
 p2 <- SIcomp_fig_fun(SI1 = SI4mix, SI2 = SI4mix_con, data = data)
 p2_empiric_nu <- SIcomp_fig_fun(SI1 = SI4mix, SI2 = SI4mix_con_empiricnu, data = data)
 p2nu0 <-  SIcomp_fig_fun(SI1 = SI4mix, SI2 = SI4mix_con_nu0, data = data)
  p4 <- SIcomp_fig_fun(SI1 = SI4mix_recall, SI2 = SI4mix_recall_con, data = data)
  p5 <- SIcomp_fig_fun(SI1 = SI3recall, SI2 = SI3recall_con, data = data)
  p6 <- SIcomp_fig_fun(SI1 = SI4recall, SI2 = SI4recall_con, data = data)

 library(cowplot)

cowplot::save_plot(filename = "figures/s3mix.jpeg", p1, base_height = 5, base_asp = 1)
 cowplot::save_plot(filename = "figures/s4mix.jpeg", p2, base_height = 5, base_asp = 1)
 cowplot::save_plot(filename = "figures/s4mix_empiric_nu.jpeg", p2_empiric_nu, base_height = 5, base_asp = 1)
 cowplot::save_plot(filename = "figures/s4mix_nu0.jpeg", p2nu0, base_height = 5, base_asp = 1)
cowplot::save_plot(filename = "figures/s4mixrecall.jpeg", p4, base_height = 5, base_asp = 1)
cowplot::save_plot(filename = "figures/s4recall.jpeg", p6, base_height = 5, base_asp = 1)
cowplot::save_plot(filename = "figures/s3recall.jpeg", p5, base_height = 5, base_asp = 1)

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


