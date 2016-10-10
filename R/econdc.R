econdc <- function(m, use_kadlec=T) {


  #compute proportion-matrix from input matrix  
  P <- pmatrix(m)
  
  #----------------------------------------
  # Test c for A conditional on B1
  
  # get data
  h1 <- P[1,1] / (P[1,1] + P[1,2])
  fa1 <- P[2,1] / (P[2,1] + P[2,2])
  h2 <- P[1,3] / (P[1,3] + P[1,4])
  fa2 <- P[2,3] / (P[2,3] + P[2,4])
  
  # run test depending on whether we use Kadlec's definition of c
  # vs. Macmillan and Creelman's definition
  if (!use_kadlec) {
    A <- testc(h1, fa1, h2, fa2, 
               sum(m[1,1:2]), sum(m[2,1:2]), sum(m[1,3:4]), sum(m[2,3:4]))
  } else {
    A <- testc2(fa1, fa2, sum(m[2,1:2]), sum(m[2,3:4]))
  }
  
  # store string with decision about test
  if (A$p_value < 0.05) {
    pass<- 'NO'
  } else {
    pass<- 'YES'
  }
  
  #store results in dataframe
  results <- data.frame(Test="c_A conditional on B1",c_hit=A$c[1], c_miss=A$c[2], z=A$z,p_value=A$p_value, Pass=pass, stringsAsFactors=FALSE)
  
  #----------------------------------------
  # Test c for A conditional on B2
  
  # get data
  h1 <- P[3,3] / (P[3,3] + P[3,4])
  fa1 <- P[4,3] / (P[4,3] + P[4,4])
  h2 <- P[3,1] / (P[3,1] + P[3,2])
  fa2 <- P[4,1] / (P[4,1] + P[4,2])
  
  # run test depending on whether we use Kadlec's definition of c
  # vs. Macmillan and Creelman's definition
  if (!use_kadlec) {
    A <- testc(h1, fa1, h2, fa2, sum(m[3,3:4]), sum(m[4,3:4]), sum(m[3,1:2]), sum(m[4,1:2]))
  } else {
    A <- testc2(fa1, fa2, sum(m[4,3:4]), sum(m[4,1:2]))
  }
  
  # store string with decision about test
  if (A$p_value < 0.05) {
    pass <- 'NO'
  } else {
    pass <- 'YES'
  }
  
  #add results to dataframe
  results<-rbind(results,c("c_A conditional on B2", c_hit=A$c[1], c_miss=A$c[2], z=A$z, p_value=A$p_value, Pass=pass))
  
  #----------------------------------------
  # Test c for B conditional on A1  
  
  # get data
  h1 <- P[1,1] / (P[1,1] + P[1,3])
  fa1 <- P[3,1] / (P[3,1] + P[3,3])
  h2 <- P[1,2] / (P[1,2] + P[1,4])
  fa2 <- P[3,2] / (P[3,2] + P[3,4])
  
  # run test depending on whether we use Kadlec's definition of c
  # vs. Macmillan and Creelman's definition
  if (!use_kadlec) {
    A <- testc(h1, fa1, h2, fa2, sum(m[1,c(1,3)]), sum(m[3,c(1,3)]), sum(m[1,c(2,4)]), sum(m[3,c(2,4)]))
  } else {
    A <- testc2(fa1, fa2, sum(m[3,c(1,3)]), sum(m[3,c(2,4)]))
  }

  # store string with decision about test
  if (A$p_value < 0.05) {
    pass <- 'NO'
  } else {
    pass <- 'YES'
  }
  
  # add results to dataframe
  results <- rbind(results, 
                   c("c_B conditional on A1", c_hit=A$c[1], c_miss=A$c[2], z=A$z, p_value=A$p_value, Pass=pass) )

  #----------------------------------------
  # Test c for B conditional on A2  

  # get data
  h1 <- P[2,2] / (P[2,2] + P[2,4])
  fa1 <- P[4,2] / (P[4,2] + P[4,4])
  h2 <- P[2,1] / (P[2,1] + P[2,3])
  fa2 <- P[4,1] / (P[4,1] + P[4,3])
  
  # run test depending on whether we use Kadlec's definition of c
  # vs. Macmillan and Creelman's definition
  if (!use_kadlec) {
    A <- testc(h1, fa1, h2, fa2, sum(m[2,c(2,4)]), sum(m[4,c(2,4)]), sum(m[2,c(1,3)]), sum(m[4,c(1,3)]))
  } else {
    A <- testc2(fa1, fa2, sum(m[4,c(2,4)]), sum(m[4,c(1,3)]))
  }
  
  # store string with decision about test
  if (A$p_value < 0.05) {
    pass <- 'NO'
  } else {
    pass <- 'YES'
  }
  
  #add results to dataframe
  results <- rbind(results, c("c_B conditional on A2", c_hit=A$c[1], c_miss=A$c[2], z=A$z, p_value=A$p_value, Pass=pass) )
  
  #output
  return(results)
}