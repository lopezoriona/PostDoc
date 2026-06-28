
# Load these auxiliary functions

circular_quantile_lsv <- function(cts, tau = 0.5) {
  
  # Median direction
  
  x_circ <- circular::circular(cts, type = "angles",
                               units = "radians", modulo = "2pi")
  circ_median <- as.numeric(circular::median.circular(x_circ))
  
  # Projection on the median
  
  theta_m_vec <- c(cos(circ_median), sin(circ_median))
  X <- cbind(cos(cts), sin(cts))
  proj <- as.numeric(X %*% theta_m_vec)
  
  # Empirical projection quantile in [-1,1]
  
  c_tau_hat <- as.numeric(stats::quantile(proj, probs = tau,
                                          type = 7, names = FALSE))
  c_tau_hat <- max(min(c_tau_hat, 1), -1)  # numerical guard
  
  # Angular distance from median
  
  delta_hat <- acos(c_tau_hat)  # in [0, pi]
  
  # Both boundary directions on the circle
  
  q_plus <- circ_median + delta_hat
  q_minus <- circ_median - delta_hat
  
  # Wrap to [0, 2*pi)
  
  q_plus <- (q_plus  + 2*pi) %% (2*pi)
  q_minus <- (q_minus + 2*pi) %% (2*pi)
  
  # Return both candidate quantiles
  
  return(c(q_plus = q_plus, q_minus = q_minus, circ_median = circ_median))
  
}


in_arc <- function(X, center, radius) {
  
  d <- abs((X - center + pi) %% (2 * pi) - pi)  
  d <= radius
  
}

