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


## Returns TOST for
## (a) parameters with maximum posterior likelihood
## (b) at the mean of the posterior distribution
## (c) a distribution of distributions of TOST ; 1
## distribution for each parameter in the posterior
## distributions of a, b, c and tmax (sampled jointly).
## tab1 is the output from fitted_params function
estimated_TOST_nf <- function(tab1, taus = seq(-20, 40, 0.1),
                              n = 1e5, fit) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- rnf(
    n = n, taus = taus, a = best["a"], b = best["b"],
    c = best["c"], tmax = best["tmax"]
  )

  mu_params <- tab1$mean
  names(mu_params) <- rownames(tab1)
  TOST_meanpars <- rnf(
    n = n,
    taus = taus,
    a = mu_params["a"],
    b = mu_params["b"],
    c = mu_params["c"],
    tmax = mu_params["tmax"]
  )

  TOST_post <- matrix(nrow = n, ncol = length(fit$a))
  for (i in 1:(length(fit$a))) {
    TOST_post[, i] <- rnf(
      n = n, taus = taus, a = fit$a[i], b = fit$b[i],
      c = fit$c[i], tmax = fit$tmax[i]
    )
  }

  TOST_post <- as.data.frame(TOST_post)

  list(
    TOST_bestpars = TOST_bestpars,
    TOST_meanpars = TOST_meanpars,
    TOST_post = TOST_post
  )
}



# sampling from empirical nu distribution instead
expected_SI_empiricnu <- function(SI, data, mixture, recall, isol, tab1, n = 1e5, tmin = -20) {
  nu <- sample(data$nu, size = n, replace = TRUE)
  SI$nu <- nu
  # data_nu <- data%>%filter(nu>0)
  # fit.gamma <- fitdist(data_nu$nu, distr = "gamma", method = "mle")
  # nu <- rgamma(n = n, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])

  # with isolation bias
  if (isol == TRUE) {
    SI <- SI %>% dplyr::filter(TOST < nu)
    SI <- sample_n(SI, n, replace = TRUE)
  }

  # with recall bias
  if (recall == TRUE) {
    recall_par <- tab1["recall", "best"]
  } else {
    recall_par <- 0
  }

  precall <- exp(-recall_par * abs(SI$SI - SI$nu))

  with_recall <- sample(
    SI$SI, n,
    replace = TRUE, prob = precall
  )

  # with invalid SIs
  if (mixture == TRUE) {
    pinvalid <- tab1["pinvalid", "best"]
  } else {
    pinvalid <- 0
  }

  toss <- runif(n)

  valid <- which(toss > pinvalid)

  invalid_si <- runif(n, tmin, 40)

  mixed <- c(
    SI$SI[valid], invalid_si[!valid]
  )

  return(pobs_SI = mixed)
}

# plot estimated SI and data
SI_figure <- function(SI, data) {
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
      label = paste(
        " mean:\n", round(mean(SI), 1), " days",
        sep = ""
      ),
      hjust = -0.1) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "bottom")
  print(p)
}

# plot to compare "true" SI and expected SI alongside the data
plot_both_SIs <- function(SI1, SI2, data) {

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
    annotate(
      geom = "text", x = mean(SI1), y = 0.085, color = "red",
      label = paste(" mean:\n", round(mean(SI1), 1), " days", sep = ""),
      hjust = -0.1
    ) +
    annotate(
      geom = "text", x = mean(SI1), y = 0.095, color = "blue",
      label = paste(" mean:\n", round(mean(SI2), 1), " days", sep = ""),
      hjust = -0.1
    ) +
    theme(legend.title = element_blank()) +
    guides(
      fill = guide_legend(override.aes = list(alpha = c(0.8, 0.3, 0.3)))
    ) +
    theme(legend.position = "bottom")


  p
}




# sampled TOST distribution --> summary statistics
TOST_summary <- function(sample) {
  mean <- mean(sample)
  sd <- sd(sample)
  ppresymp <- 100 * sum(sample < 0) / length(sample)
  p10days <- 100 * sum(sample < 10) / length(sample)
  lower <- quantile(sample, 0.01)
  upper <- quantile(sample, 0.99)

  c(
    mean_inf = mean,
    sd_inf = sd,
    presymp_inf = ppresymp,
    after10_inf = p10days,
    l_quantile_inf = lower,
    u_quantile_inf = upper
  )
}

# put summary stats from SI and TOST into a table (for a single model)
tost_si_summary <- function(tost_samples,si_samples, digits = 2) {
  foo <- as.data.frame(apply(tost_samples$TOST_post, 2, TOST_summary))
  CrI_2.5 <- apply(foo, 1, quantile, probs = 0.025)
  CrI_97.5 <- apply(foo, 1, quantile, probs = 0.975)

  tab2 <- data.frame(
    best_pars = TOST_summary(tost_samples$TOST_bestpars),
    mean_pars = TOST_summary(tost_samples$TOST_meanpars),
    CrI_2.5 = CrI_2.5,
    CrI_97.5 = CrI_97.5
  )
  mean_si <- c(
    mean(si_samples$SI_bestpars$SI),
    mean(si_samples$SI_meanpars$SI),
    quantile(
      colMeans(si_samples$SI_post), probs = c(0.025, 0.975),
      na.rm = TRUE
    )
  )
  median_si <- c(
    median(si_samples$SI_bestpars$SI),
    median(si_samples$SI_meanpars$SI),
    quantile(
      colMedians(as.matrix(si_samples$SI_post)),
      probs = c(0.025, 0.975), na.rm = TRUE
    )
  )
  sd_si <- c(
    sd(si_samples$SI_bestpars$SI),
    sd(si_samples$SI_meanpars$SI),
    quantile(
      colSds(as.matrix(si_samples$SI_post)),
      probs = c(0.025, 0.975), na.rm = TRUE
    )
  )

  tab_2 <- rbind(
    tab2, mean_si = mean_si, median_si = median_si,
    sd_si = sd_si
  )

  round(tab_2, digits)
}

# plot TOST
TOST_figure <- function(TOST) {
  cutoff <- 0
  hist.y <- density(TOST, from = -20, to = 40) %$%
    data.frame(x = x, y = y) %>%
    mutate(area = x >= cutoff)
  ## also print area left of cutoff
  presymp <- scales::percent(sum(TOST < 0) / length(TOST), 0.1)
  label <- glue::glue("Pre-symptomatic\n{presymp}")
  the.plot <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) +
    geom_ribbon() +
    geom_line(aes(y = y)) +
    geom_vline(xintercept = cutoff, color = "red") +
    annotate(
      geom = "text", x = cutoff, y = max(hist.y$y) * 0.95,
      color = "red", label = "symptom\nonset", hjust = -0.1
    ) +
    annotate(
      geom = "text", x = cutoff - 19,
      y = max(hist.y$y) * 0.9,
      color = "red",
      label = label, hjust = -0.1
    ) +
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
expected_SI <- function(SI, data, mixture, recall, isol, tab1, n = 1e5, tmin = -20) {
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
SI_figure <- function(SI, data) {
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
plot_both_SIs <- function(SI1, SI2, data) {
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

  tab2 <- tost_si_summary(samples)

  p1 <- TOST_figure(samples$TOST_bestpars)

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

  p2 <- plot_both_SIs(SI1 = SI, SI2 = SI_conditional, data = data)

  return(list(
    table1 = tab1,
    table2 = tab2,
    p1 = p1,
    p2 = p2
  ))
}
