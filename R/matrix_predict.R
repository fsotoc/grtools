matrix_predict<-function(means,covmat,b,c){
  # computes a predicted matrix of confusion probabilities
  # means is a 4x2 matrix in which each row vector represents the mean of one perceptual distribution
  # covmat is a list of 4 2x2 covariance matrices, one for each perceptual distribution
  # b is a 2x2 matrix with one row for each vector of slope parameters
  # c is a column vector (2x1 matrix) with constants from decision bounds
  
  
  z_vector = c(-2.5758,   -2.1900,   -1.9673,   -1.8157,   -1.6978,   -1.5998,   -1.5153,   -1.4404,   -1.3729 ,  -1.3112 ,  -1.2540 ,  -1.2008,   -1.1507,
               -1.1034,   -1.0584,   -1.0154,   -0.9743,   -0.9348 ,  -0.8966,   -0.8598 ,  -0.8240,   -0.7893 ,  -0.7555  , -0.7226 ,  -0.6904,   -0.6589,
               -0.6281,   -0.5978,   -0.5681,   -0.5389,   -0.5101 ,  -0.4816,   -0.4538 ,  -0.4262,   -0.3989 ,  -0.3719,   -0.3452 ,  -0.3167,   -0.2924,
               -0.2663,   -0.2404,   -0.2147,   -0.1691,   -0.1637 ,  -0.1383,   -0.1130 ,  -0.0879,   -0.0627 ,  -0.0376,   -0.0125 ,   0.0125,    0.0376,
               0.0627,    0.0879,   0.1130 ,   0.1383  ,  0.1637 ,   0.1691  ,  0.2147   , 0.2404  ,  0.2663   , 0.2924  ,  0.3167   , 0.3452  ,  0.3719,
               0.3989,    0.4262,    0.4538,    0.4816 ,   0.5101,    0.5389 ,   0.5681  ,  0.5978 ,   0.6281  ,  0.6589 ,   0.6904  ,  0.7226 ,   0.7555,
               0.7893,    0.8240,    0.8598,    0.8966 ,   0.9348,    0.9743 ,   1.0154  ,  1.0584 ,   1.1034  ,  1.1507 ,   1.2008  ,  1.2540 ,   1.3112,
               1.3729,    1.4404,    1.5153,    1.5998 ,   1.6978,    1.8157 ,   1.9673  ,  2.1900 ,   2.5758);
  
  #compute matrix of p(r_j|s_i)
  p_matrix <- matrix(0, nrow=4, ncol=4, byrow=T)
  z <- matrix(nrow=2, ncol=1, byrow=TRUE)
  for (i in 1:4){
    #compute p(r_j|s_i)
    
    #compute P_i from Cholesky factorization
    P <-matrix(0,2,2)
    P[1,1]<-sqrt(covmat[[i]][1,1])
    P[2,1]<-covmat[[i]][2,1]/P[1,1]
    P[2,2]<-sqrt(covmat[[i]][2,2]-P[2,1]^2)
    
    # get constant values for discriminant function
    scal1<-b[1,]%*%P
    scal2<-b[2,]%*%P
    const1<-b[1,]%*%means[i,]+c[1,]
    const2<-b[2,]%*%means[i,]+c[2,]
    
    
    p<-matrix(0,1,4)+.000001
    for(k in 1:100){
      z[1]<-z_vector[k]
      for(l in 1:100){
        z[2]<-z_vector[l]
        #apply both discriminant functions to z-vector
        h_1<-scal1%*%z+const1
        h_2<-scal2%*%z+const2
        if (h_1 < 0){
          if (h_2< 0){
            p[1]<-p[1]+.0001
          } else {
            p[3]<-p[3]+.0001}
        } else {
          if (h_2 < 0){
            p[2]<-p[2]+.0001
          } else {
            p[4]<-p[4]+.0001
          }
        }
      }
    }
    p_matrix[i,] <- p/sum(p)
  }
  
  return(p_matrix)
  
}