simulate_4a_mix <- function(params_inc, params_inf, params_iso,
                            offset, max_shed, pinvalid, nsim,
                            alpha_invalid, beta_invalid, min_si, max_si,
                            p_recall) {
  
  nsim_local <- 10*nsim
  
  valid_si <- better_simulate_si(
    params_inc, params_inf, params_iso,
    min_inf = offset, max_inf = max_shed, nsim = nsim_local)
  
  valid_si <- valid_si%>%
    dplyr::filter(t_1<nu)%>%
    dplyr::select(si, nu)
  valid_si$group <- "valid"
  
  valid_si <- valid_si[sample(nrow(valid_si), nsim), ]
  
  invalid_si <- data.frame(
    si = rbeta(nsim, shape1 = alpha_invalid, shape2 = beta_invalid),
    nu = rgamma(nsim, shape = params_iso$shape, scale = params_iso$scale)
  )
  f <- map_into_interval(0, 1, min_si, max_si)
  invalid_si$si <- f(invalid_si$si)
  invalid_si$group <- "invalid"
  
  simulated_si <- valid_si[0, ]
  
  for (idx in seq_len(nsim)) {
    toss <- runif(1)
    if (toss <= pinvalid) {
      simulated_si <- rbind(
        simulated_si, invalid_si[sample(nrow(invalid_si), 1), ]
      )
    } else {
      simulated_si <- rbind(
        simulated_si, valid_si[sample(nrow(valid_si), 1), ]
      )
    }
  }
  
  # sample using probability of recall
  simulated_si$toss2 <- runif(nsim, min = 0, max = 1)
  
  simulated_si_recalled <- simulated_si%>%
    dplyr::mutate(recall = exp(-abs(si-nu)*p_recall))%>%
    dplyr::filter(toss2 < recall)
  
  list(
    valid_si = valid_si,
    invalid_si = invalid_si,
    simulated_si = simulated_si,
    simulated_si_recalled = simulated_si_recalled
  )
}

