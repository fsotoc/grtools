testp <- function(pp1, pp2, n1, n2) {
  
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