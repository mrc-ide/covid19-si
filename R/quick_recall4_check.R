recall4_1 <- readRDS("4a_recall_sim_1.rds")
recall4_2 <- readRDS("4a_recall_sim_2.rds")
recall4_3 <- readRDS("4a_recall_sim_3.rds")
recall4_1_nonorm <- readRDS("stanfits/4a_recall_sim_no_norm1.rds")
nsim <- 3
recall41 <- rstan::extract(recall4_1_nonorm)
#rstan::summary(recall41)

idx <- which.max(recall41[["lp__"]])
best_params <- list(lp__ = recall41[["lp__"]][idx],
                    alpha = recall41[["alpha1"]][idx],
                    beta = recall41[["beta1"]][idx],
                    recall = recall41[["recall"]][idx])
offset <- -3
max_shed <- 21
mean <- ((max_shed - offset) *
  beta_shape1shape22muvar(shape1 = best_params[["alpha"]],
        shape2 = best_params[["beta"]])[["mu"]])+offset

recov_inf <- mu_sd_posterior_distr(recall41, 21, -3)

mean_best3 <- recov_inf[which((recov_inf$var == "best") & (recov_inf$param == "mu")),]$val
sd_best3 <- recov_inf[which((recov_inf$var == "best") & (recov_inf$param == "sd")),]$val
recall_best3 <- best_params[["recall"]]

true_mean <- c(2,8,2)
true_sd <- c(1,2,1)
true_recall <- c(0.4, 0.4, 0.4)

df <- data.frame(mean = c(mean_best1,mean_best2,
                             mean_best3),
                    sd = c(sd_best1, sd_best2, sd_best3),
                    recall = c(recall_best1, recall_best2, recall_best3), 
                 sim = 1:nsim)
dflong <- reshape2::melt(df, id.var = "sim")
truedf <- data.frame(true_mean, true_sd, true_recall, sim = 1:nsim)
truedflong <- reshape2::melt(truedf, id.var = "sim", value.name = "true_value")
dflong$true_value <- truedflong$true_value

ggplot(data = dflong)+
  geom_point(aes(x = sim, y = value), shape = 16)+
  geom_point(aes(x = sim, y = true_value), shape = 4)+
  facet_wrap(~variable, nrow = 3, scales = "free")
