#' @importFrom Rcpp evalCpp
#' @useDynLib grtools
grt_wind_nll_r3x2 <- function(data,params){
  
   N <- length(data)
   L <- 0
   
   # means are the same for everybody
   means<-matrix(params[1:8],4,2,byrow=TRUE)
   
   covmat<-list()
   
   # compute log-likelihood for each subject and add to total value
   for (sub in 1:N){
     
     # get individual attention parameters
     kap = params[20+(sub-1)*8+1]
     lam1 = params[20+(sub-1)*8+2]
     lam2 = 1-lam1
          
     # scale the covariance matrices for this subject
     covmat[[1]] <- matrix(c(params[9]/(kap*lam1), params[11]*sqrt((params[9]/(kap*lam1))*(params[10]/(kap*lam2))),
                             params[11]*sqrt((params[9]/(kap*lam1))*(params[10]/(kap*lam2))), params[10]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[2]] <- matrix(c(params[12]/(kap*lam1), params[14]*sqrt((params[12]/(kap*lam1))*(params[13]/(kap*lam2))),
                             params[14]*sqrt((params[12]/(kap*lam1))*(params[13]/(kap*lam2))), params[13]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[3]] <- matrix(c(params[15]/(kap*lam1), params[17]*sqrt((params[15]/(kap*lam1))*(params[16]/(kap*lam2))),
                             params[17]*sqrt((params[15]/(kap*lam1))*(params[16]/(kap*lam2))), params[16]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[4]] <- matrix(c(params[18]/(kap*lam1), params[20]*sqrt((params[18]/(kap*lam1))*(params[19]/(kap*lam2))),
                             params[20]*sqrt((params[18]/(kap*lam1))*(params[19]/(kap*lam2))), params[19]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     # get decision bound parameters for this subject
     b <- matrix(c(1, params[20+(sub-1)*8+3],
                   1, params[20+(sub-1)*8+5],
                   params[20+(sub-1)*8+7], 1),
                  3,2, byrow=TRUE)
     c <- matrix(c(params[20+(sub-1)*8+4], 
                   params[20+(sub-1)*8+6],
                   params[20+(sub-1)*8+8]),3,1)
     
     L <- L + matrixloglikC_r3x2(data[[sub]],means,covmat,b,c)
     
   }

   L <- -L
   return(L)
   
}

