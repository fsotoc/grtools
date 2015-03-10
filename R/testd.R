testd <- function(h1, fa1, h2, fa2, n1_n, n1_s, n2_n, n2_s){
  
  # compute d1 and d2
  d1 <- qnorm(h1) - qnorm(fa1)
  d2 <- qnorm(h2) - qnorm(fa2)

  # compute variance of d1 and d2
  var_d1 <- ( (fa1 * (1 - fa1) ) / (n1_n * dnorm( qnorm(fa1) )^2) ) +
    ( (h1 * (1 - h1)) / (n1_s * dnorm( qnorm(h1) )^2) )
  var_d2 <- ( (fa2 * (1 - fa2) ) / (n2_n * dnorm( qnorm(fa2) )^2) ) +
    ( (h2 * (1 - h2)) / (n2_s * dnorm( qnorm(h2) )^2) )

  #compute z
  z <- (d1 - d2) / (sqrt(var_d1 + var_d2))
  
  #compute p-value, two-tailed
  p_val <- 2 * (pnorm(-abs(z)))
  
  #output
  return(list(d=(c(d1,d2)), var=(c(var_d1,var_d2)), z=z, p_value=p_val))
}