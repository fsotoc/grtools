#' @importFrom Rcpp evalCpp
#' @useDynLib grtools
grt_wind_nll <- function(data,params){
  
   N <- length(data)
   L <- 0
   
   # means are the same for everybody
   means<-matrix(c(0,0,params[1:6]),4,2,byrow=TRUE)
   
   covmat<-list()
   
   # compute log-likelihood for each subject and add to total value
   for (sub in 1:N){
     
     # get individual attention parameters
     kap = params[16+(sub-1)*6+1]
     lam1 = params[16+(sub-1)*6+2]
     lam2 = 1-lam1
          
     # scale the covariance matrices for this subject
     covmat[[1]] <- matrix(c(1/(kap*lam1), params[7]*(1/sqrt(lam1*lam2*kap^2)),
                             params[7]*(1/sqrt(lam1*lam2*kap^2)), 1/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[2]] <- matrix(c(params[8]/(kap*lam1), params[10]*sqrt((params[8]/(kap*lam1))*(params[9]/(kap*lam2))),
                             params[10]*sqrt((params[8]/(kap*lam1))*(params[9]/(kap*lam2))), params[9]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[3]] <- matrix(c(params[11]/(kap*lam1), params[13]*sqrt((params[11]/(kap*lam1))*(params[12]/(kap*lam2))),
                             params[13]*sqrt((params[11]/(kap*lam1))*(params[12]/(kap*lam2))), params[12]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[4]] <- matrix(c(params[14]/(kap*lam1), params[16]*sqrt((params[14]/(kap*lam1))*(params[15]/(kap*lam2))),
                             params[16]*sqrt((params[14]/(kap*lam1))*(params[15]/(kap*lam2))), params[15]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     # get decision bound parameters for this subject
     b <- matrix(c(1, params[16+(sub-1)*6+3], params[16+(sub-1)*6+5], 1),
                  2,2,byrow=TRUE)
     c <- matrix(c(params[16+(sub-1)*6+4], params[16+(sub-1)*6+6]),2,1)
     
     L <- L + matrixloglikC(data[[sub]],means,covmat,b,c)
     
   }

   L <- -L
   return(L)
   
}
