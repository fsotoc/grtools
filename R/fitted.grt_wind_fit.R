#' @export
fitted.grt_wind_fit <- function(model, out_format="vector"){
  # input (model) is a list obtained from grt_wind_fit
  # out_format determines the format of the output; "vector" returns
  # a vector of predicted response probabilities; "list" returns
  # a list of confusion probability matrices (one per participant)
  
   N <- model$N
   predicted <- c()
   
   # means are the same for everybody
   means<-matrix(c(0,0,model$par[1:6]),4,2,byrow=TRUE)
   
   covmat<-list()
   
   # compute log-likelihood for each subject and add to total value
   for (sub in 1:N){
     
     # get individual attention parameters
     kap = model$par[16+(sub-1)*6+1]
     lam1 = model$par[16+(sub-1)*6+2]
     lam2 = 1-lam1
          
     # scale the covariance matrices for this subject
     covmat[[1]] <- matrix(c(1/(kap*lam1), model$par[7]*(1/sqrt(lam1*lam2*kap^2)),
                             model$par[7]*(1/sqrt(lam1*lam2*kap^2)), 1/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[2]] <- matrix(c(model$par[8]/(kap*lam1), model$par[10]*sqrt((model$par[8]/(kap*lam1))*(model$par[9]/(kap*lam2))),
                             model$par[10]*sqrt((model$par[8]/(kap*lam1))*(model$par[9]/(kap*lam2))), model$par[9]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[3]] <- matrix(c(model$par[11]/(kap*lam1), model$par[13]*sqrt((model$par[11]/(kap*lam1))*(model$par[12]/(kap*lam2))),
                             model$par[13]*sqrt((model$par[11]/(kap*lam1))*(model$par[12]/(kap*lam2))), model$par[12]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[4]] <- matrix(c(model$par[14]/(kap*lam1), model$par[16]*sqrt((model$par[14]/(kap*lam1))*(model$par[15]/(kap*lam2))),
                             model$par[16]*sqrt((model$par[14]/(kap*lam1))*(model$par[15]/(kap*lam2))), model$par[15]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     # get decision bound parameters for this subject
     b <- matrix(c(1, model$par[16+(sub-1)*6+3], model$par[16+(sub-1)*6+5], 1),
                  2,2,byrow=TRUE)
     c <- matrix(c(model$par[16+(sub-1)*6+4], model$par[16+(sub-1)*6+6]),2,1)
     
     if (out_format=="vector"){
       predicted <- c(predicted, as.vector(matrix_predict(means,covmat,b,c)))
     } else if (out_format=="list"){
       predicted <- c(predicted, list(matrix_predict(means,covmat,b,c)))
     }
   }

   return(predicted)
   
}
