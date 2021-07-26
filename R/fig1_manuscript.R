# FIGURE 1C

# read in simulated dataset
sim_data <- readRDS("data/sim_data_manuscript.rds")

# s3 for figure 1
fit3_sim <- readRDS("data/illustration_s3.rds")
tab1_s3 <- fitted_params(fit3_sim)
tost_s3 <- estimated_TOST_nf(tab1 = tab1_s3, 
                             fit = rstan::extract(fit3_sim))
si_s3 <- unconditional_si(inf_samples = tost_s3$TOST_bestpars)
si_s3_df <- data.frame(SI = si_s3$si, tag = "baseline")

# s4 for figure 1
fit4_sim <- readRDS("data/illustration_s4.rds")
tab1_s4 <- fitted_params(fit4_sim)
tost_s4 <- estimated_TOST_nf(tab1 = tab1_s4, 
                             fit = rstan::extract(fit4_sim))
si_s4 <- unconditional_si(inf_samples = tost_s4$TOST_bestpars)
si_s4_df <- data.frame(SI = si_s4$si, tag = "isolation")

# reshaping data for ggplot
sim_data_m <- rbind(sim_data$SI_full, sim_data$SI_observed) %>%
  dplyr::select(SI, tag)
data_gg <- rbind(sim_data_m, si_s3_df, si_s4_df)
data_gg$tag[which(data_gg$tag == "obs")] <- "observed"
data_gg$tag[which(data_gg$tag == "full")] <- "true"


p1 <- ggplot()+
  geom_histogram(data = data_gg%>%filter(tag == "observed" | tag == "true"), 
                 aes(x = SI, y = ..density.., fill = tag), position = "identity", color = NA, binwidth = 1, alpha = 0.6)+
  scale_fill_manual(values = c("black", "#999999"), name = "Simulated SI data")+  
  geom_density(data = data_gg%>%filter(tag == "baseline" | tag == "isolation"), 
               aes(x = SI, y = ..density.., lty = tag), color = "#D55E00", fill = NA, size = 1)+
  scale_linetype_manual(values = c(1,2), name = "Model")+
  theme_minimal()+
  xlab("Serial Interval (days)")+
  theme(legend.position = "bottom")+
  xlim(NA, 50)
  
cowplot::save_plot("1c_v2.tiff", p1, base_width = 6.5, base_height = 2.75)
