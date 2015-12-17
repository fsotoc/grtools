testc <- function(h1, fa1, h2, fa2, n1_n, n1_s, n2_n, n2_s) {
  
  # replace 0 and 1 by close approximations to avoid Inf values
  props = c(h1, fa1, h2, fa2)
  if (any(props==0)){warning("Proportions in the data equal to 0; replaced by 1e-10")}
  props[props==0] <- 1e-10
  if (any(props==1)){warning("Proportions in the data equal to 1; replaced by 1-(1e-10)")}
  props[props==1] <- 1-1e-10
  h1 <- props[1]
  fa1 <- props[2]
  h2 <- props[3]
  fa2 <- props[4]
  
  # compute c1 and c2
  c1 <- (qnorm(h1) + qnorm(fa1)) / 2
  c2 <- (qnorm(h2) + qnorm(fa2)) / 2
  
  # compute variance of c1 and c2
  var_c1 <- ( (fa1*(1-fa1)) / (n1_n*dnorm(qnorm(fa1))^2) ) + 
    ( (h1*(1-h1)) / (n1_s*dnorm(qnorm(h1))^2) ) / 4
  var_c2 <- ( (fa2*(1-fa2)) / (n2_n*dnorm(qnorm(fa2))^2) ) + 
    ( (h2*(1-h2)) / (n2_s*dnorm(qnorm(h2))^2) ) / 4

  # compute z
  z <- (c1 - c2) / (sqrt(var_c1 + var_c2))

  #compute p-value, two-tailed
  p_val <- 2*(pnorm(-abs(z)))
  
  #output
  return(list(c=(c(c1,c2)), z=z, p_value=p_val))
}