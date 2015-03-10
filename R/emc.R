emc<-function(m, use_kadlec=T){
  #compute proportion-matrix from input matrix  
  P<-pmatrix(m)
  
  #compute c for A across the two levels of B
  h1 <- P[1,1] + P[1,3]
  fa1 <- P[2,1] + P[2,3]
  h2 <- P[3,1] + P[3,3]
  fa2 <- P[4,1] + P[4,3]

  # compute c and test for significance
  if (!use_kadlec) {
    A <- testc(h1, fa1, h2, fa2, sum(m[1,]), sum(m[2,]), sum(m[3,]), sum(m[4,]))
  } else {
    A <- testc2(fa1, fa2, sum(m[2,]), sum(m[4,]))
  }

  # check whether test was significant
  if (A$p_value < 0.05) {
    Pass <- 'NO'
  } else {
    Pass <- 'YES'
  }
  
  #store results in dataframe
  results <- data.frame(Test="c_A across the two levels of B",
                        c1=A$c[1], c2=A$c[2], z=A$z, 
                        p_value=A$p_value, Pass=Pass, stringsAsFactors=FALSE)
  
  #compute c for B across the two levels of A
  h1 <- P[1,1] + P[1,2]
  fa1 <- P[3,1] + P[3,2]
  h2 <- P[2,1] + P[2,2]
  fa2 <- P[4,1] + P[4,2]

  # compute c and test for significance
  if (!use_kadlec) {
    A <- testc(h1, fa1, h2, fa2, sum(m[2,]), sum(m[1,]), sum(m[4,]), sum(m[3,]))
  } else {
    A <- testc2(fa1, fa2, sum(m[2,]), sum(m[4,]))
  }
  
  # check whether test was significant
  if (A$p_value < 0.05) {
    Pass<- 'NO'
  } else {
    Pass<- 'YES'
  }
  
  #add results to dataframe
  results <- rbind(results,c("c_B across the two levels of A",
                           c1=A$c[1], c2=A$c[2], z=A$z,
                           p_value=A$p_value, Pass=Pass))
  
  #output
  return(results)
}
