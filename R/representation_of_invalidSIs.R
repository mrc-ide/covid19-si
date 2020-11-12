simulate_3a_mix <- function(mean_inc, sd_inc, shape1, shape2, pinvalid, nsim, offset,
                            alpha_invalid, beta_invalid, min_si, max_si) {
  valid_si <- (simulate_si(
    mean_inc, sd_inc, shape1, shape2,
    max_shed, NULL, NULL, nsim = nsim
  )$si) - offset
  
  max_si <- max(valid_si)
  
  invalid_si <- (max_si - min_si)*
    rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid)
  
  invalid_si <- invalid_si + min_si
  toss <- runif(nsim)
  
  sim_invalid <- data.frame(si = invalid_si[toss <= pinvalid], group = "invalid")
  sim_valid <- data.frame(si = valid_si[toss > pinvalid], group = "valid")
  
  simulated_si <- rbind(sim_invalid, sim_valid) 
  
  list(
    valid_si = valid_si,
    invalid_si = invalid_si,
    simulated_si = simulated_si
  )
}

si_sim_test <- simulate_3a_mix(mean_inc_og, sd_inc_og, shape1_max, 
                shape2_max, p_invalid_max, nsim = 100000, offset = offset,
                alpha_invalid, beta_invalid, min_si = min(data_c$si),
                max_si = max(data_c$si))

si_sim_test$simulated_si$group <- factor(si_sim_test$simulated_si$group, levels = c("valid",
                                                                                    "invalid"))


psi <- ggplot(data = si_sim_test$simulated_si, aes(fill = factor(group))) +
  geom_density(aes(si_sim_test$simulated_si$si),
               alpha = 0.3,
               position = "stack"
  ) 

psi <- ggplot(data = si_sim_test$simulated_si, aes(fill = factor(group))) +
  geom_histogram(aes(si_sim_test$simulated_si$si),
               alpha = 0.3,
               position = "stack",
               binwidth = 1
  ) 



+
  geom_histogram(
    data = data_c, aes(si, y = ..density.., fill = "blue"),
    alpha = 0.3,
    binwidth = 1
  ) +
  geom_vline(
    xintercept = mean(data_c$si), col = "blue", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(si_sim_test$simulated_si$si), col = "red", linetype = "dashed"
  ) +
  scale_fill_identity(
    guide = "legend",
    labels = c("Data", "Posterior (valid SIs"),
    breaks = c("blue", "red")
  ) +
  theme_minimal() +
  xlab("Serial Interval") +
  theme(legend.title = element_blank())
