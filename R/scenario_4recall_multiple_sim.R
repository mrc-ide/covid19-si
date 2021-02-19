## Set up params grid
prefix <- "4a_recall_sim_no_norm"

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = "inc_par2",
  params_iso = c("iso_par1", "iso_par2"),
  params_offset = "offset3", # just try with offset = -3 initially
  params_pinvalid = "pinvalid1", #just try with pinvalid = 0 initially
  params_recall = c("recall1", "recall2"),
  stringsAsFactors = FALSE
)

index <- 1:nrow(param_grid)
param_grid <- param_grid[index, ]

params_inf_all <- pmap(
  list(
    params_inf = param_grid$params_inf,
    params_offset = param_grid$params_offset
  ),
  function(params_inf, params_offset) {
    out <- params[[params_inf]]
    ## The whole shifting in simulation will shift mu,
    ## so pass larger my to simulate function/
    offset <- params[[params_offset]]
    beta_muvar2shape1shape2(
      (out$mean_inf - offset) / (max_shed - offset),
      out$sd_inf^2 / (max_shed - offset)^2
    )
  }
)

params_inc_all <- map(
  param_grid$params_inc,
  function(params_inc) {
    out <- params[[params_inc]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_iso_all <- map(
  param_grid$params_iso,
  function(params_iso) {
    out <- params[[params_iso]]
    epitrix::gamma_mucv2shapescale(
      mu = out[[1]], cv = out[[2]] / out[[1]]
    )
  }
)

params_offset_all <- map(
  param_grid$params_offset,
  function(params_offset) params[[params_offset]]
)

params_pinvalid_all <- map(
  param_grid$params_pinvalid,
  function(params_pinvalid) params[[params_pinvalid]]
)

params_recall_all <- map(
  param_grid$params_recall,
  function(params_recall) params[[params_recall]]
)

simulated <- pmap(
  list(
    params_inc = params_inc_all,
    params_inf = params_inf_all,
    params_iso = params_iso_all,
    params_offset = params_offset_all,
    params_pinvalid = params_pinvalid_all,
    params_recall = params_recall_all
  ),
  function(params_inc, params_inf, params_iso, params_offset,
           params_pinvalid, params_recall) {
    min_si <- params_offset
    simulate_4a_mix2(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_post_filter, alpha_invalid, beta_invalid,
      min_si, max_si, params_recall
    )
  }
)

sampled <- map(
  simulated, function(df) {
    ## might have NA from when pinvalid = 0
    out <- df[["simulated_si_recalled"]]
    out <- out[complete.cases(out), ]
    #idx <- sample(nsim_pre_filter, nsim_post_filter, replace = TRUE) - I sampled in my sim function
    out
  }
)

#si_vec <- seq(-3, 40, by = 1)

fits <- pmap(
  list(
    params_inc = params_inc_all,
    params_offset = params_offset_all,
    sim_data = sampled,
    index = index
  ),
  function(params_inc, params_offset, sim_data, index) {
    ## Rounding now to check things
    sim_data$si <- round(sim_data$si)
    fit_4a_recall <- stan(
      file = here::here("stan-models/full_model_no_norm.stan"),
      data = list(
        N = length(sim_data$si),
        si = sim_data$si,
        nu = sim_data$nu,
        max_shed = max_shed,
        offset1 = params_offset,
        max_si = max_si,
        min_si = params_offset - 0.001,  
        alpha2 = params_inc[["shape"]],
        beta2 = 1 / params_inc[["scale"]],
        width = width
        #M = length(si_vec),
        #y_vec = si_vec
      ),
      chains = 2,
      iter = 2000,
      seed = 42,
      verbose = TRUE
    )
    ## Save it now in case R crashes
    saveRDS(fit_4a_recall, glue::glue("stanfits/{prefix}{index}.rds"))
    fit_4a_recall
  }
)

## infiles <- glue::glue("stanfits/{prefix}{index}.rds")
## fits <- map(infiles, readRDS)

process_fits <- pmap_dfr(
  list(
    fit = fits,
    true_vals = param_grid$params_inf,
    params_inc = param_grid$params_inc,
    params_iso = param_grid$params_iso,
    offset = params_offset_all,
    pinvalid = params_pinvalid_all),
  function(fit, true_vals, params_inc, params_iso, offset, pinvalid) {
    
    true_vals <- params[[true_vals]]
    incubation <- params[[params_inc]]
    isolation <- params[[params_iso]]
    
    true_val_df <- data.frame(
      true_val = c(true_vals$mean_inf, true_vals$sd_inf, pinvalid),
      param = c("mu", "sd", "pinvalid")
    )
    samples <- rstan::extract(fit)
    pinvalid_posterior <- quantile_as_df(samples[["pinvalid"]])
    pinvalid_posterior$param <- "pinvalid"
    out <- mu_sd_posterior_distr(samples, max_shed, offset)
    out <- rbind(out, pinvalid_posterior) %>%
      tidyr::spread(var, val)
    
    out <- left_join(out, true_val_df)
    out$mean_inc <- incubation$mean_inc
    out$sd_inc <- incubation$sd_inc
    out$mean_iso <- isolation$mean_iso
    out$sd_iso <- isolation$sd_iso
    out
  }, .id = "sim"
)

outfile <- glue::glue("data/{prefix}params_posterior_distr.rds")
saveRDS(process_fits, outfile)


## process_fits <- readRDS(outfile)

process_fits <- mutate_if(process_fits, is.numeric, ~ signif(., 3))
##process_fits$sim <- as.integer(process_fits$sim)

x <- split(process_fits, process_fits$param)

names(x) <- c("Mean (Infectious Profile)", "pinvalid",
              "SD (infectious profile)")

iwalk(
  x,
  function(df, param) {
    df$mean_iso <- as.factor(df$mean_iso)
    df <- arrange(df, mean_iso)
    df$sim <- factor(df$sim, levels = df$sim, ordered = TRUE)
    
    p <- ggplot(df) +
      geom_linerange(
        aes(
          x = sim, ymin = `25%`, ymax = `75%`, col = mean_iso
        )
      ) +
      geom_point(aes(sim, `50%`, col = mean_iso)) +
      geom_point(aes(sim, true_val), shape = 4) +
      facet_grid(mean_inc ~ true_val, scales = "free") +
      expand_limits(y = 0) +
      theme_minimal() +
      theme(
        legend.position = "top", axis.text.x = element_blank(),
        axis.title.x = element_blank()
      ) +
      ylab(param)
    
    ggsave(glue::glue("figures/{prefix}{param}.png"))
  }
)

##################################################
### SI posterior distribution
##################################################
best_params <- map(
  fits,
  function(fit) {
    out <- rstan::extract(fit)
    idx <- which.max(out[["lp__"]])
    map(out, ~ .[idx])
  }
)

params_inf_post <- map(
  best_params, ~ list(shape1 = .$alpha1, shape2 = .$beta1)
)

params_pinvalid_post <- map(best_params, ~ .$pinvalid)

post_simulated <- pmap(
  list(
    params_inc = params_inc_all,
    params_inf = params_inf_post,
    params_iso = params_iso_all,
    params_offset = params_offset_all,
    params_pinvalid = params_pinvalid_post
  ),
  function(params_inc, params_inf, params_iso, params_offset,
           params_pinvalid) {
    min_si <- params_offset
    simulate_3a_mix(
      params_inc, params_inf, params_iso, params_offset, max_shed,
      params_pinvalid, nsim_pre_filter, alpha_invalid, beta_invalid,
      min_si, max_si
    )
  }
)

post_sampled <- map(
  post_simulated, function(df) {
    ## might have NA from when pinvalid = 0
    out <- df[["simulated_si"]]
    out <- out[complete.cases(out), ]
    idx <- sample(nsim_pre_filter, nsim_post_filter, replace = TRUE)
    out[idx, ]
  }
)

post_si_plots <- map2(
  sampled, post_sampled, function(x, y) {
    ggplot() +
      geom_density(aes(round(x$si), fill = "red"), alpha = 0.3, col = NA) +
      geom_density(aes(round(y$si), fill = "blue"), alpha = 0.3, col = NA) +
      scale_fill_identity(
        breaks = c("red", "blue"),
        labels = c("SI used for fitting", "SI Posterior distributio"),
        guide = "legend"
      ) +
      xlab("Serial Interval") +
      theme_minimal() +
      theme(legend.position = "top", legend.title = element_blank())
    
  }
)

outfiles <- glue::glue("figures/{prefix}{index}_si.png")

walk2(
  post_si_plots, outfiles, function(p, filename) ggsave(filename, p)
)


post_si_qntls <- map2_dfr(
  sampled, post_sampled, function(x, y) {
    before <- quantile_as_df(x$si)
    before$category <- "used"
    after <- quantile_as_df(y$si)
    after$category <- "posterior"
    rbind(before, after)
  }, .id = "sim")


post_si_qntls <- tidyr::spread(post_si_qntls, var, val)
##post_si_qntls$sim <- as.integer(post_si_qntls$sim)
post_si_qntls <- arrange(post_si_qntls, `50%`)
post_si_qntls$sim <- factor(
  post_si_qntls$sim,
  levels = post_si_qntls$sim[! duplicated(post_si_qntls$sim)],
  ordered = TRUE
)

ggplot(post_si_qntls) +
  geom_point(
    aes(sim, `50%`, col = category),
    position = position_dodge(width = 0.5)) +
  geom_linerange(
    aes(x = sim, ymin = `25%`, ymax = `75%`, col = category),
    position = position_dodge(width = 0.5)) +
  ylab("Serial Interval") +
  theme_minimal() +
  theme(
    legend.position = "top", axis.text.x = element_blank(),
    axis.title.x = element_blank(), legend.title = element_blank()
  )

outfile <- glue::glue("figures/{prefix}serial_interval.png")

ggsave(outfile)
