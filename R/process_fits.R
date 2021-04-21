############################
## process & present fits ##
############################
#functions
      dnf <- function(t, a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
          numerator <- a + b
          growing <- a * exp(b * (t - tmax))
          failing <- b * exp(-a * (t - tmax))
          (numerator / (growing + failing))^c
      }
      
      rnf <- function(n, taus = seq(-20, 40, 0.1), a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
          p <- map_dbl(taus, function(t) dnf(t, a, b, c, tmax))
          sampled <- sample(taus, n, replace = TRUE, prob = p)
      }
        
        #a. extract fitted parameters
      table1_fun <- function(fit){
          tab1 <- as.data.frame(rstan::summary(fit)$summary)
          fit <- rstan::extract(fit)
          
          best_idx <- which.max(fit[["lp__"]])
          best <- unlist(map(fit, function(x) x[[best_idx]]))
          tab1 <- add_column(tab1, best = best, .after = 1)
          return(tab1)
      }
        
        #b. fitted parameters --> sampled distributions
      sample_dist_fun <- function(tab1, taus = seq(-20, 40, 0.1), n = 1e4, fit,
                                    shape_inc = 1.675513, scale_inc = 3.043843){
          # TOST
          best <- tab1$best
          names(best) <- rownames(tab1)
          TOST_bestpars <- rnf(n = n,
                               taus = taus, 
                               a = best["a"],
                               b = best["b"], 
                               c = best["c"], 
                               tmax = best["tmax"])
          
          mean <- tab1$mean
          names(mean) <- rownames(tab1)
          TOST_meanpars <- rnf(n = n,
                               taus = taus,
                               a = mean["a"],
                               b = mean["b"],
                               c = mean["c"],
                               tmax = mean["tmax"])
          
          TOST_post <- matrix(ncol = length(fit$a),
                              nrow = n)
          for(i in 1:(length(fit$a))){
            TOST_post[,i] <- rnf(n = n,
                                 taus = taus,
                                 a = fit$a[i],
                                 b = fit$b[i],
                                 c = fit$c[i],
                                 tmax = fit$tmax[i])         
          }
          
          TOST_post <- as.data.frame(TOST_post)
          
          #SI
          inc_mat <- as.data.frame(matrix(rgamma(shape = shape_inc, 
                                                 scale = scale_inc, 
                                                 n = n*length(fit$a)),
                                          nrow = dim(TOST_post)[1], 
                                          ncol = dim(TOST_post)[2]))
          SI_best <- inc_mat[,1] + TOST_bestpars
          SI_mean <- inc_mat[,1] + TOST_meanpars
          SI_post <- TOST_post + inc_mat
          
          return(list(TOST_bestpars = TOST_bestpars,
                      TOST_meanpars = TOST_meanpars,
                      TOST_post = TOST_post,
                      SI_bestpars = SI_best,
                      SI_meanpars = SI_mean,
                      SI_post = SI_post))
          
      }
        
        #c. sampled distributions --> summary statistics
      TOST_summary_func <- function(sample){
          mean <- mean(sample)
          sd <- sd(sample)
          ppresymp <- 100*sum(sample<0)/length(sample)
          p10days <- 100*sum(sample<10)/length(sample)
          lower <- quantile(sample, 0.01)
          upper <- quantile(sample, 0.99)
          
          return(c(mean_inf = mean,
                   sd_inf = sd,
                   presymp_inf = ppresymp,
                   after10_inf = p10days,
                   l_quantile_inf = lower,
                   u_quantile_inf = upper))
      }
      
      table2_fun <- function(samples){
        foo <- as.data.frame(apply(samples$TOST_post, 2, TOST_summary_func))
        CrI_2.5 <- apply(foo, 1, quantile, probs = 0.025)
        CrI_97.5 <- apply(foo, 1, quantile, probs = 0.975)
        
        tab2 <- data.frame(best_pars = TOST_summary_func(samples$TOST_bestpars),
                           mean_pars = TOST_summary_func(samples$TOST_meanpars),
                           CrI_2.5 = CrI_2.5,
                           CrI_97.5 = CrI_97.5)
        mean_si <- c(mean(samples$SI_bestpars), mean(samples$SI_meanpars), 
                     quantile(colMeans(samples$SI_post), probs = c(0.025, 0.975)))
        median_si <- c(median(samples$SI_bestpars), median(samples$SI_meanpars), 
                       quantile(colMedians(as.matrix(samples$SI_post)), probs = c(0.025, 0.975)))
        sd_si <- c(sd(samples$SI_bestpars), sd(samples$SI_meanpars), 
                   quantile(colSds(as.matrix(samples$SI_post)), probs = c(0.025, 0.975)))
        tab_2 <- rbind(tab2, mean_si = mean_si, median_si = median_si, sd_si = sd_si)
        
        return(tab_2)
      }
        
      TOST_fig_fun <- function(TOST3, TOST4){

        tost_df <- data.frame(id = 1:(length(TOST3)),
                              No = TOST3,
                              Yes = TOST4)
        tost_df <- reshape2::melt(tost_df, id.vars = "id")
        tost_df <- tost_df%>%
          mutate(presymp = value<0)
        ggplot(data = tost_df)+
          geom_density(aes(value, group = variable, fill = presymp), size = 1)+
          xlab("TOST (days)")+
          theme_minimal()+
          labs(colour = "corrected\nisolation bias")+
          geom_vline(xintercept = 0, lty = 2)
      }
      
      SI_fig_fun <- function(SI3, SI4, data){
        
      }
        
  # wrapper function
        

  #1. read in rds files

  fit3 <- readRDS("stanfits/release/scenario3a_mixture_nf.rds")
  fit3_recall <- readRDS("stanfits/release/scenario3arecall_mixture_nf.rds")
  fit4 <- readRDS("stanfits/release/scenario4a_mixture_nf.rds")
  fit4_recall <- readRDS("stanfits/release/scenario4arecall_mixture_nf.rds")

  
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
  
  p <- TOST_fig_fun(TOST3, TOST4)
  
  SI3 <- samples_s3$SI_bestpars
  SI4 <- samples_s4$SI_bestpars
  data <- readRDS("data/cowling_data_clean.rds")
  data <- data%>%
    mutate(si = as.numeric(si))%>%
    dplyr::rename(nu = onset_first_iso)
  