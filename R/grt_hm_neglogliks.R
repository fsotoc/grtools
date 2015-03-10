# note here how the likelihood function requires both w and data as inputs
# otherwise it wouldn't be general (it would work with only one data table)

#### logL_PI_PS_DS ####
negloglik_mod1<-function(w,data){ 
  # w is just a vector with 4 values: two means (w(1) and w(2) below) and the position
  # of two decision bounds (w(3) and w(4))
    
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[3,2]<-w[2]
  means[4,1]<-w[1]
  means[4,2]<-w[2]
 
  
  # get covariance matrices (fix variances to 1)
  covmat<-list()
  covmat[[1]]<-diag(2)
  covmat[[2]]<-diag(2)
  covmat[[3]]<-diag(2)
  covmat[[4]]<-diag(2)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[3],w[4]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_PI_PS_A_DS #### 
negloglik_mod2<-function(w,data){ 
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[2,2]<-w[2]
  means[3,2]<-w[3]
  means[4,1]<-w[1]
  means[4,2]<-w[4]
  
  
  # get covariance matrices (fix variances to 1)
  covmat<-list()
  covmat[[1]]<-diag(2)
  covmat[[2]]<-diag(2)
  covmat[[3]]<-diag(2)
  covmat[[4]]<-diag(2)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[5],w[6]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_PI_PS_B_DS #### 
negloglik_mod3<-function(w,data){
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[3,1]<-w[2]
  means[3,2]<-w[3]
  means[4,1]<-w[4]
  means[4,2]<-w[3]
  
  # get covariance matrices (fix variances to 1)
  covmat<-list()
  covmat[[1]]<-diag(2)
  covmat[[2]]<-diag(2)
  covmat[[3]]<-diag(2)
  covmat[[4]]<-diag(2)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[5],w[6]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_1RHO_PS_DS #### 
negloglik_mod4<-function(w,data){
  
  # get means
  
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[3,2]<-w[2]
  means[4,1]<-w[1]
  means[4,2]<-w[2]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  vals<-c(1,w[3],w[3],1)
  covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[4],w[5]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_PS_A_1RHO_DS #### 
negloglik_mod5<-function(w,data){
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[2,2]<-w[2]
  means[3,2]<-w[3]
  means[4,1]<-w[1]
  means[4,2]<-w[4]
  
  # get covariance matrices (fix variances to 1)
  covmat<-list()
  vals<-c(1,w[5],w[5],1)
  covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[6],w[7]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_PI_DS #### 
negloglik_mod6<-function(w,data){
  
  # get means
  
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[2,2]<-w[2]
  means[3,1]<-w[3]
  means[3,2]<-w[4]
  means[4,1]<-w[5]
  means[4,2]<-w[6]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  covmat[[1]]<-diag(2)
  covmat[[2]]<-diag(2)
  covmat[[3]]<-diag(2)
  covmat[[4]]<-diag(2)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[7],w[8]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}


#### logL_PS_B_1RHO_DS ####
negloglik_mod7<-function(w,data){
  
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[3,1]<-w[2]
  means[3,2]<-w[3]
  means[4,1]<-w[4]
  means[4,2]<-w[3]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  vals<-c(1,w[5],w[5],1)
  covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[6],w[7]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_PS_DS ####
negloglik_mod8<-function(w,data){
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[3,2]<-w[2]
  means[4,1]<-w[1]
  means[4,2]<-w[2]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  covmat[[1]]<-matrix(c(1,w[3],w[3],1),2,2,byrow=TRUE)
  covmat[[2]]<-matrix(c(1,w[4],w[4],1),2,2,byrow=TRUE)
  covmat[[3]]<-matrix(c(1,w[5],w[5],1),2,2,byrow=TRUE)
  covmat[[4]]<-matrix(c(1,w[6],w[6],1),2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[7],w[8]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}


#### logL_PS_A_DS ####
negloglik_mod9<-function(w,data){
  
  # get means
  
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[2,2]<-w[2]
  means[3,2]<-w[3]
  means[4,1]<-w[1]
  means[4,2]<-w[4]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  covmat[[1]]<-matrix(c(1,w[5],w[5],1),2,2,byrow=TRUE)
  covmat[[2]]<-matrix(c(1,w[6],w[6],1),2,2,byrow=TRUE)
  covmat[[3]]<-matrix(c(1,w[7],w[7],1),2,2,byrow=TRUE)
  covmat[[4]]<-matrix(c(1,w[8],w[8],1),2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[9],w[10]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_1RHO_DS ####
negloglik_mod10<-function(w,data){
  
  # get means
  
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[2,2]<-w[2]
  means[3,1]<-w[3]
  means[3,2]<-w[4]
  means[4,1]<-w[5]
  means[4,2]<-w[6]
  
  # get covariance matrices (fix variances to 1)
  covmat<-list()
  vals<-c(1,w[7],w[7],1)
  covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
  covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[8],w[9]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}

#### logL_PS_B_DS #### 
negloglik_mod11<-function(w,data){
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[3,1]<-w[2]
  means[3,2]<-w[3]
  means[4,1]<-w[4]
  means[4,2]<-w[3]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  covmat[[1]]<-matrix(c(1,w[5],w[5],1),2,2,byrow=TRUE)
  covmat[[2]]<-matrix(c(1,w[6],w[6],1),2,2,byrow=TRUE)
  covmat[[3]]<-matrix(c(1,w[7],w[7],1),2,2,byrow=TRUE)
  covmat[[4]]<-matrix(c(1,w[8],w[8],1),2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[9],w[10]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}


#### logL_DS #### 
negloglik_mod12<-function(w,data){
  
  # get means
  means<-matrix(0,4,2,byrow=TRUE)
  means[2,1]<-w[1]
  means[2,2]<-w[2]
  means[3,1]<-w[3]
  means[3,2]<-w[4]
  means[4,1]<-w[5]
  means[4,2]<-w[6]
  
  # get covariance matrices (fix variances to 1)
  
  covmat<-list()
  covmat[[1]]<-matrix(c(1,w[7],w[7],1),2,2,byrow=TRUE)
  covmat[[2]]<-matrix(c(1,w[8],w[8],1),2,2,byrow=TRUE)
  covmat[[3]]<-matrix(c(1,w[9],w[9],1),2,2,byrow=TRUE)
  covmat[[4]]<-matrix(c(1,w[10],w[10],1),2,2,byrow=TRUE)
  
  # get decision bound parameters
  b<-diag(2)
  c<-matrix(c(w[11],w[12]),2,1)
  
  L<-matrixloglikC(data,means,covmat,b,c)
  return(-L)
}
