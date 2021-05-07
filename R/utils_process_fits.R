############################
## process & present fits ##
############################

# functions
# density function for NF distribution
dnf <- function(t, a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
  numerator <- a + b
  growing <- a * exp(b * (t - tmax))
  failing <- b * exp(-a * (t - tmax))
  (numerator / (growing + failing))^c
}

# random sample from NF distribution
rnf <- function(n, taus = seq(-20, 40, 0.1), a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
  p <- map_dbl(taus, function(t) dnf(t, a, b, c, tmax))
  sample(taus, n, replace = TRUE, prob = p)
}


      # predicted observed SIs (under assumed biases)
      # currently assumes recall and isolation biases only affect valid SIs
      expected_SI_fun <- function(SI, data, mixture, recall, isol, tab1, n = 1e5, tmin = -20){

        # fit.gamma <- fitdist(data$nu - tmin, distr = "gamma", method = "mle")
        # nu_shifted <- rgamma(n = n, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])
        # nu <- nu_shifted + tmin
        # SI$nu <- nu
        data_nu <- data%>%filter(nu>0)
        fit.gamma <- fitdist(data_nu$nu, distr = "gamma", method = "mle")
        nu <- rgamma(n = n, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])
        SI$nu <- nu

        # with isolation bias
        if(isol == TRUE) {
          SI <- SI%>%dplyr::filter(TOST<nu)
          SI <- sample_n(SI, n, replace = TRUE)
        }

        # with recall bias
        if(recall == TRUE) {
          recall_par <- tab1["recall", "best"]
        } else recall_par <- 0

        precall <- exp(-recall_par * abs(SI$SI - SI$nu))

        with_recall <- sample(
          SI$SI, n, replace = TRUE, prob = precall
        )

        # with invalid SIs
        if(mixture == TRUE) {
          pinvalid <- tab1["pinvalid", "best"]
        } else pinvalid <- 0

        toss <- runif(n)

        valid <- which(toss > pinvalid)

        invalid_si <- runif(n, tmin, 40)

        mixed <- c(
          SI$SI[valid], invalid_si[!valid]
        )

       return(pobs_SI = mixed)
      }

      # sampling from empirical nu distribution instead
      expected_SI_fun_empiricnu <- function(SI, data, mixture, recall, isol, tab1, n = 1e5, tmin = -20){

        nu <- sample(data$nu, size = n, replace = TRUE)
        SI$nu <- nu
        # data_nu <- data%>%filter(nu>0)
        # fit.gamma <- fitdist(data_nu$nu, distr = "gamma", method = "mle")
        # nu <- rgamma(n = n, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])

        # with isolation bias
        if(isol == TRUE) {
          SI <- SI%>%dplyr::filter(TOST<nu)
          SI <- sample_n(SI, n, replace = TRUE)
        }

        # with recall bias
        if(recall == TRUE) {
          recall_par <- tab1["recall", "best"]
        } else recall_par <- 0

        precall <- exp(-recall_par * abs(SI$SI - SI$nu))

        with_recall <- sample(
          SI$SI, n, replace = TRUE, prob = precall
        )

        # with invalid SIs
        if(mixture == TRUE) {
          pinvalid <- tab1["pinvalid", "best"]
        } else pinvalid <- 0

        toss <- runif(n)

        valid <- which(toss > pinvalid)

        invalid_si <- runif(n, tmin, 40)

        mixed <- c(
          SI$SI[valid], invalid_si[!valid]
        )

        return(pobs_SI = mixed)
      }

      # plot estimated SI and data
      SI_fig_fun <- function(SI, data){
        p <-  ggplot() +
          geom_histogram(
            data = data, aes(si, y = ..density.., fill = "gray77"),
            alpha = 0.8,
            binwidth = 1
          ) +
          geom_density(aes(SI, fill = "red"),
                       alpha = 0.3, colour = NA
          ) +
          geom_vline(
            xintercept = mean(data$si), col = "gray77", linetype = "dashed"
          ) +
          geom_vline(
            xintercept = mean(SI), col = "red", linetype = "dashed"
          ) +
          scale_fill_identity(
            guide = "legend",
            labels = c("Data", "Posterior"),
            breaks = c("gray77", "red")
          ) +
          theme_minimal() +
          xlab("Serial Interval (days)") +
          annotate(geom = 'text', x = mean(SI), y = 0.095, color = 'red', label = paste(" mean:\n",round(mean(SI), 1), " days", sep = ""), hjust = -0.1)+
          theme(legend.title = element_blank())+
          theme(legend.position="bottom")
        print(p)
      }

      # plot to compare "true" SI and expected SI alongside the data
      SIcomp_fig_fun <- function(SI1, SI2, data){
        p <-  ggplot() +
          geom_histogram(
            data = data, aes(si, y = ..density.., fill = "gray77"),
            alpha = 0.8,
            binwidth = 1
          ) +
          geom_density(aes(SI1, fill = "red"),
                       alpha = 0.3, colour = NA
          ) +
          geom_density(aes(SI2, fill = "blue"),
                       alpha = 0.3, colour = NA
          ) +
          geom_vline(
            xintercept = mean(data$si), col = "gray77", linetype = "dashed"
          ) +
          geom_vline(
            xintercept = mean(SI1), col = "red", linetype = "dashed"
          ) +
          geom_vline(
            xintercept = mean(SI2), col = "blue", linetype = "dashed"
          ) +
          scale_fill_identity(
            guide = "legend",
            labels = c("Data", "Posterior - True", "Posterior - Expected"),
            breaks = c("gray77", "red", "blue")
          ) +
          theme_minimal() +
          xlab("Serial Interval (days)") +
          annotate(geom = 'text', x = mean(SI1), y = 0.085, color = 'red', label = paste(" mean:\n",round(mean(SI1), 1), " days", sep = ""), hjust = -0.1)+
          annotate(geom = 'text', x = mean(SI1), y = 0.095, color = 'blue', label = paste(" mean:\n",round(mean(SI2), 1), " days", sep = ""), hjust = -0.1)+
          theme(legend.title = element_blank())+
          guides(fill = guide_legend(override.aes= list(alpha = c(0.8, 0.3,0.3))))+
          theme(legend.position="bottom")
        print(p)
      }


# extract fitted parameters
table1_fun <- function(fit) {
  tab1 <- as.data.frame(rstan::summary(fit)$summary)
  fit <- rstan::extract(fit)

  best_idx <- which.max(fit[["lp__"]])
  best <- unlist(map(fit, function(x) x[[best_idx]]))
  tab1 <- add_column(tab1, best = best, .after = 1)

  tab1 <- round(tab1, 2)
  return(tab1)
}

# fitted parameters --> sampled TOST and SI distributions
sample_dist_fun <- function(tab1, taus = seq(-20, 40, 0.1), n = 1e5, fit,
                            shape_inc = params_real$inc_par2[["shape"]],
                            scale_inc = params_real$inc_par2[["scale"]]) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- rnf(
    n = n,
    taus = taus,
    a = best["a"],
    b = best["b"],
    c = best["c"],
    tmax = best["tmax"]
  )

  mean <- tab1$mean
  names(mean) <- rownames(tab1)
  TOST_meanpars <- rnf(
    n = n,
    taus = taus,
    a = mean["a"],
    b = mean["b"],
    c = mean["c"],
    tmax = mean["tmax"]
  )

  TOST_post <- matrix(
    ncol = length(fit$a),
    nrow = n
  )
  for (i in 1:(length(fit$a))) {
    TOST_post[, i] <- rnf(
      n = n,
      taus = taus,
      a = fit$a[i],
      b = fit$b[i],
      c = fit$c[i],
      tmax = fit$tmax[i]
    )
  }

  TOST_post <- as.data.frame(TOST_post)

  # SI
  inc_mat <- as.data.frame(matrix(rgamma(
    shape = shape_inc,
    scale = scale_inc,
    n = n * length(fit$a)
  ),
  nrow = dim(TOST_post)[1],
  ncol = dim(TOST_post)[2]
  ))
  SI_best <- inc_mat[, 1] + TOST_bestpars
  SI_mean <- inc_mat[, 1] + TOST_meanpars
  SI_post <- TOST_post + inc_mat

  return(list(
    TOST_bestpars = TOST_bestpars,
    TOST_meanpars = TOST_meanpars,
    TOST_post = TOST_post,
    SI_bestpars = data.frame(SI = SI_best, TOST = TOST_bestpars),
    SI_meanpars = data.frame(SI = SI_mean, TOST = TOST_meanpars),
    SI_post = SI_post
  ))
}

# sampled TOST distribution --> summary statistics
TOST_summary_func <- function(sample) {
  mean <- mean(sample)
  sd <- sd(sample)
  ppresymp <- 100 * sum(sample < 0) / length(sample)
  p10days <- 100 * sum(sample < 10) / length(sample)
  lower <- quantile(sample, 0.01)
  upper <- quantile(sample, 0.99)

  return(c(
    mean_inf = mean,
    sd_inf = sd,
    presymp_inf = ppresymp,
    after10_inf = p10days,
    l_quantile_inf = lower,
    u_quantile_inf = upper
  ))
}

# put summary stats from SI and TOST into a table (for a single model)
table2_fun <- function(samples) {
  foo <- as.data.frame(apply(samples$TOST_post, 2, TOST_summary_func))
  CrI_2.5 <- apply(foo, 1, quantile, probs = 0.025)
  CrI_97.5 <- apply(foo, 1, quantile, probs = 0.975)

  tab2 <- data.frame(
    best_pars = TOST_summary_func(samples$TOST_bestpars),
    mean_pars = TOST_summary_func(samples$TOST_meanpars),
    CrI_2.5 = CrI_2.5,
    CrI_97.5 = CrI_97.5
  )
  mean_si <- c(
    mean(samples$SI_bestpars$SI), mean(samples$SI_meanpars$SI),
    quantile(colMeans(samples$SI_post), probs = c(0.025, 0.975))
  )
  median_si <- c(
    median(samples$SI_bestpars$SI), median(samples$SI_meanpars$SI),
    quantile(colMedians(as.matrix(samples$SI_post)), probs = c(0.025, 0.975))
  )
  sd_si <- c(
    sd(samples$SI_bestpars$SI), sd(samples$SI_meanpars$SI),
    quantile(colSds(as.matrix(samples$SI_post)), probs = c(0.025, 0.975))
  )
  tab_2 <- rbind(tab2, mean_si = mean_si, median_si = median_si, sd_si = sd_si)

  tab_2 <- round(tab_2, 2)

  return(tab_2)
}

# plot TOST
TOST_fig_fun <- function(TOST) {
  cutoff <- 0
  hist.y <- density(TOST, from = -20, to = 40) %$%
    data.frame(x = x, y = y) %>%
    mutate(area = x >= cutoff)

  the.plot <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) +
    geom_ribbon() +
    geom_line(aes(y = y)) +
    geom_vline(xintercept = cutoff, color = "red") +
    annotate(geom = "text", x = cutoff, y = max(hist.y$y) * 0.95, color = "red", label = "symptom\nonset", hjust = -0.1) +
    theme_minimal() +
    xlab("TOST (days)") +
    ylab("density") +
    theme(legend.position = "none") +
    scale_fill_manual(values = alpha("gray", c(.3, 0)))
  print(the.plot)
}

# plot to compare TOSTS
TOSTcomp_fig_fun <- function(TOST3, TOST4) {
  tost_df <- data.frame(
    id = 1:(length(TOST3)),
    No = TOST3,
    Yes = TOST4
  )
  tost_df <- reshape2::melt(tost_df, id.vars = "id")
  tost_df <- tost_df %>%
    mutate(presymp = value < 0)
  p <- ggplot(data = tost_df) +
    geom_density(aes(value, group = variable, fill = presymp), size = 1) +
    xlab("TOST (days)") +
    theme_minimal() +
    labs(colour = "corrected\nisolation bias") +
    geom_vline(xintercept = 0, lty = 2)

  print(p)
}

# predicted observed SIs (under assumed biases)
# currently assumes recall and isolation biases only affect valid SIs
expected_SI_fun <- function(SI, data, mixture, recall, isol, tab1, n = 1e5, tmin = -20) {
  fit.gamma <- fitdist(data$nu - tmin, distr = "gamma", method = "mle")
  nu_shifted <- rgamma(n = n, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])
  nu <- nu_shifted + tmin

  # data_nu <- data%>%filter(nu>0)
  # fit.gamma <- fitdist(data_nu$nu, distr = "gamma", method = "mle")
  # nu <- rgamma(n = n, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])

  # with isolation bias
  if (isTRUE(isol)) {
    SI$nu <- nu
    SI <- SI %>% dplyr::filter(TOST < nu)
    SI <- sample_n(SI, n * 0.1, replace = TRUE)
  }

  # with recall bias
  if (isTRUE(recall)) {
    recall_par <- tab1["recall", "best"]
  } else {
    recall_par <- 0
  }

  precall <- exp(-recall_par * abs(SI$SI - SI$nu))

  with_recall <- sample(
    SI$SI, n * 0.1,
    replace = TRUE, prob = precall
  )

  # with invalid SIs
  if (isTRUE(mixture)) {
    pinvalid <- tab1["pinvalid", "best"]
  } else {
    pinvalid <- 0
  }

  toss <- runif(n * 0.1)

  valid <- which(toss > pinvalid)

  invalid_si <- runif(n * 0.1, tmin, 40)

  mixed <- c(
    SI$SI[valid], invalid_si[!valid]
  )

  return(pobs_SI = mixed)
}

# plot estimated SI and data
SI_fig_fun <- function(SI, data) {
  p <- ggplot() +
    geom_histogram(
      data = data, aes(si, y = ..density.., fill = "gray77"),
      alpha = 0.8,
      binwidth = 1
    ) +
    geom_density(aes(SI, fill = "red"),
      alpha = 0.3, colour = NA
    ) +
    geom_vline(
      xintercept = mean(data$si), col = "gray77", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(SI), col = "red", linetype = "dashed"
    ) +
    scale_fill_identity(
      guide = "legend",
      labels = c("Data", "Posterior"),
      breaks = c("gray77", "red")
    ) +
    theme_minimal() +
    xlab("Serial Interval (days)") +
    annotate(
      geom = "text", x = mean(SI), y = 0.095, color = "red",
      label = paste(" mean:\n", round(mean(SI), 1), " days", sep = ""),
      hjust = -0.1
    ) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "bottom")
  print(p)
}

# plot to compare "true" SI and expected SI alongside the data
SIcomp_fig_fun <- function(SI1, SI2, data) {
  p <- ggplot() +
    geom_histogram(
      data = data, aes(si, y = ..density.., fill = "gray77"),
      alpha = 0.8,
      binwidth = 1
    ) +
    geom_density(aes(SI1, fill = "red"),
      alpha = 0.3, colour = NA
    ) +
    geom_density(aes(SI2, fill = "blue"),
      alpha = 0.3, colour = NA
    ) +
    geom_vline(
      xintercept = mean(data$si), col = "gray77", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(SI1), col = "red", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(SI2), col = "blue", linetype = "dashed"
    ) +
    scale_fill_identity(
      guide = "legend",
      labels = c("Data", "Posterior - True", "Posterior - Expected"),
      breaks = c("gray77", "red", "blue")
    ) +
    theme_minimal() +
    xlab("Serial Interval (days)") +
    annotate(geom = "text", x = mean(SI1), y = 0.085, color = "red", label = paste(" mean:\n", round(mean(SI1), 1), " days", sep = ""), hjust = -0.1) +
    annotate(geom = "text", x = mean(SI1), y = 0.095, color = "blue", label = paste(" mean:\n", round(mean(SI2), 1), " days", sep = ""), hjust = -0.1) +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(override.aes = list(alpha = c(0.8, 0.3, 0.3)))) +
    theme(legend.position = "bottom")
  print(p)
}

wrapper_single_model <- function(stanfit, data, mixture = T, recall = F, isol = F, n = 1e4) {
  tab1 <- table1_fun(fit = stanfit)

  samples <- sample_dist_fun(tab1, fit = rstan::extract(stanfit))

  tab2 <- table2_fun(samples)

  p1 <- TOST_fig_fun(samples$TOST_bestpars)

  SI <- samples$SI_bestpars$SI

  SI_conditional <- pred_observed_SI_fun(
    SI = samples$SI_bestpars,
    data = data,
    mixture = mixture,
    recall = recall,
    isol = isol,
    tab1 = tab1,
    n = n
  )

  p2 <- SIcomp_fig_fun(SI1 = SI, SI2 = SI_conditional, data = data)

  return(list(
    table1 = tab1,
    table2 = tab2,
    p1 = p1,
    p2 = p2
  ))
}
