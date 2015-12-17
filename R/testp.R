testp <- function(pp1, pp2, n1, n2) {
  
  # replace 0 and 1 by close approximations to avoid Inf values
  props = c(pp1, pp2)
  if (any(props==0)){warning("Proportions in the data equal to 0; replaced by 1e-10")}
  props[props==0] <- 1e-10
  if (any(props==1)){warning("Proportions in the data equal to 1; replaced by 1-(1e-10)")}
  props[props==1] <- 1-1e-10
  pp1 <- props[1]
  pp2 <- props[2]

  #compute z-score
  if (pp1==pp2){
    z <- 0
  } else {
    z<- (pp1 - pp2) / ( sqrt((pp1 * (1 - pp1) / n1) + 
                              (pp2 * (1-pp2) / n2)) )
  }
  
  #compute p-value, two-tailed
  p_val <- 2 * (pnorm(-abs(z)))
  
  #output
  return(list(z=z, p_value=p_val))
}