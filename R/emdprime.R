emdprime <- function(m){
  
  #compute proportion-matrix from input matrix  
  P <- pmatrix(m)
  
  #------------------------------------------
  # Test d' for A across the two levels of B
  
  # get data
  h1 <- P[1,1] + P[1,3]
  fa1 <- P[2,1] + P[2,3]
  h2 <- P[3,1] + P[3,3]
  fa2 <- P[4,1] + P[4,3]
  
  #compute d' and test for significant difference
  A <- testd(h1, fa1, h2, fa2, sum(m[1,]), sum(m[2,]), sum(m[3,]), sum(m[4,]))
  
  # check whether test was significant
  if (A$p_value < 0.05) {
    Pass <- 'NO'
  } else {
    Pass<- 'YES'
  }
  
  #store results in dataframe
  results <- data.frame(Test="d'_A across the two levels of B",
                        d1=A$d[1], d2=A$d[2], z=A$z, 
                        p_value=A$p_value, Pass=Pass, stringsAsFactors=FALSE)
  
  #------------------------------------------
  # Test d' for B across the two levels of A
  
  # get data
  h1 <- P[1,1] + P[1,2]
  fa1 <- P[3,1] + P[3,2]
  h2 <- P[2,1] + P[2,2]
  fa2 <- P[4,1] + P[4,2]
  
  # compute d' and test for significance
  A <- testd(h1, fa1, h2, fa2, sum(m[1,]), sum(m[3,]), sum(m[2,]), sum(m[4,]))
  
  # check whether test was significant
  if (A$p_value < 0.05) {
    Pass <- 'NO'
  } else {
    Pass <- 'YES'
  }
  
  #add results to dataframe
  results <- rbind( results, c("d'_B across the two levels of A", 
                              d1=A$d[1], d2=A$d[2], z=A$z,
                              p_value=A$p_value, Pass=Pass) )

  #output
  return(results)
}
