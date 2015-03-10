testc2 <- function(fa1, fa2, n1_n, n2_n){
  
  # compute c1 and c2
  c1 <- qnorm(fa1)
  c2 <- qnorm(fa2)
  
  # compute variance of c1 and c2
  var_c1 <- ( (fa1*(1 - fa1)) / (n1_n*dnorm(qnorm(fa1))^2) )
  var_c2 <- ( (fa2*(1 - fa2)) / (n2_n*dnorm(qnorm(fa2))^2) )

  #compute z-score
  z <- (c1 - c2) / (sqrt(var_c1 + var_c2))
  
  #compute p-value, two-tailed
  p_val <- 2*(pnorm(-abs(z)))
  
  #output
  return(list(c=(c(c1,c2)), var=(c(var_c1,var_c2)), z=z, p_value=p_val))
}