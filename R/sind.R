sind <- function(m){
  #compute proportion-matrix from input matrix  
  P <- pmatrix(m)
  
  #create matrix
  index <- matrix(c(1,1,2,1,1,2,2,2), nrow=4, ncol=2, byrow=TRUE)
  
  #create lists with the stimulus names and response names
  stimtxt <- c('A1B1', 'A2B1', 'A1B2', 'A2B2')
  resptxt <- c('a1b1', 'a2b1', 'a1b2', 'a2b2')
  
  #computing marginals, by row, eliciting [A1, A2, B1, B2]
  for (i in 1:4) {
    marg_p <- matrix(c((P[i,1] + P[i,3]),
                       (P[i,2] + P[i,4]),
                       (P[i,1] + P[i,2]),
                       (P[i,3] + P[i,4])),
                     nrow=4, ncol=2, byrow=TRUE)
    
    #computing expected probability (ep) from sampling indepedence (SI)
    for (j in 1:4){    
      ep <- marg_p[1,index[j,1]] * marg_p[2,index[j,2]]
      
      # test proportions
      A <- testp(P[i,j], ep, sum(m[1,]), sum(m[1,]))
      
      # check whether test was significant
      if (A$p_value < 0.05)   {
        Pass<- 'NO'
      } else {
        Pass<- 'YES'
      }
      
      # store results in dataframe    
      if (i==1 && j==1) {
        results <- data.frame(Stimulus=stimtxt[i], Response=resptxt[j], 
                              Expected_p=ep, Observed_p=P[i,j], 
                              z=A$z, p_value=A$p_value, Pass=Pass, 
                              stringsAsFactors=FALSE)
      } else {
        results <- rbind(results, c(Stimulus=stimtxt[i], Response=resptxt[j],
                                    Expected_p=ep, observed_p=P[i,j], 
                                    z=A$z, p_value=A$p_value, Pass=Pass))
      }  
    }
  }
  
  #output
  return(results)
}
