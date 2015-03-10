mri_test<-function(m){

  # compute proportion-matrix from input matrix  
  P <- pmatrix(m)
  
  # compute p(a_1|A_1) across the two levels of B
  p_1 <- P[1,1] + P[1,3]
  p_2 <- P[3,1] + P[3,3]
  
  # perform test on proportions
  A <- testp(p_1, p_2, sum(m[1,]) ,sum(m[3,]))
  
  # check whether the test was significant
  if (A$p_value < 0.05) {
    Pass<- 'NO'
  } else {
    Pass<- 'YES'
  }
  
  # store results in dataframe
  results <- data.frame(Test="p(a_1|A_1) across the two levels of B",
                      z=A$z, p_value=A$p_value, Pass=Pass, stringsAsFactors=FALSE)
  
  # compute p(a_2|A_2) across the two levels of B
  p_1 <- P[2,2]+P[2,4]
  p_2 <- P[4,2]+P[4,4]
  
  # perform test on proportions
  A <- testp(p_1, p_2,sum(m[2,]), sum(m[4,]))
  
  # check whether the test was significant
  if (A$p_value < 0.05) {
    Pass<- 'NO'
  } else {
    Pass<- 'YES'
  }
  
  # add results to dataframe
  results <- rbind(results,c("p(a_2|A_2) across the two levels of B",
                             A$z, A$p_value, Pass))
  
  # compute p(b_1|B_1) across the two levels of A
  p_1 <- P[1,1] + P[1,2]
  p_2 <- P[2,1] + P[2,2]
  
  # perform test on proportions
  A <- testp(p_1, p_2, sum(m[1,]), sum(m[2,]))
  
  # check whether the test was significant
  if (A$p_value < 0.05) {
    Pass<- 'NO'
  } else {
    Pass<- 'YES'
  }
  
  # add results to dataframe
  results <- rbind(results,c("p(b_1|B_1) across the two levels of A",
                             A$z, A$p_value, Pass))
  
  # compute p(b_2|B_2) across the two levels of A
  p_1 <- P[3,3] + P[3,4]
  p_2 <- P[4,3] + P[4,4]
  
  # perform test on proportions
  A <- testp(p_1, p_2, sum(m[3,]), sum(m[4,]))
  
  # check whether the test was significant
  if (A$p_value < 0.05) {
    Pass<- 'NO'
  } else {
    Pass<- 'YES'
  }
  
  # add results to dataframe
  results <- rbind(results,c("p(b_2|B_2) across the two levels of A",
                             A$z, A$p_value, Pass))
  # output
  return(results)
}
  