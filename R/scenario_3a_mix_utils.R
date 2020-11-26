simulate_3a_mix <- function(params_inc, params_inf, params_iso, 
                            offset, max_shed, pinvalid, nsim, 
                            alpha_invalid, beta_invalid, min_si, max_si) {
  
  valid_si <- better_simulate_si(
    params_inc, params_inf, params_iso, 
    min_inf = offset, max_inf = max_shed, nsim)$si
  
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
