simulate_3a_mix <- function(mean_inc, sd_inc, shape1, shape2, max_shed,
                            pinvalid, nsim, offset, alpha_invalid, 
                            beta_invalid, min_si, max_si) {
  valid_si <- (simulate_si(
    mean_inc, sd_inc, shape1, shape2,
    max_shed, NULL, NULL, nsim = nsim
  )$si) + offset
  
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
