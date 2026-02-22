#' @export
fitted_grt_wind_fit_r3x2 <- function(model, out_format="vector"){
  # input (model) is a list obtained from grt_wind_fit
  # out_format determines the format of the output; "vector" returns
  # a vector of predicted response probabilities; "list" returns
  # a list of confusion probability matrices (one per participant)
  
   N <- model$N
   predicted <- c()
   
   # means are the same for everybody
   means<-matrix(model$fullpars[1:8],4,2,byrow=TRUE)
   
   covmat<-list()
   
   # compute log-likelihood for each subject and add to total value
   for (sub in 1:N){
     
     # get individual attention parameters
     kap = model$fullpars[20+(sub-1)*8+1]
     lam1 = model$fullpars[20+(sub-1)*8+2]
     lam2 = 1-lam1
          
     # scale the covariance matrices for this subject
     covmat[[1]] <- matrix(c(model$fullpars[9]/(kap*lam1), model$fullpars[11]*sqrt((model$fullpars[9]/(kap*lam1))*(model$fullpars[10]/(kap*lam2))),
                             model$fullpars[11]*sqrt((model$fullpars[9]/(kap*lam1))*(model$fullpars[10]/(kap*lam2))), model$fullpars[10]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[2]] <- matrix(c(model$fullpars[12]/(kap*lam1), model$fullpars[14]*sqrt((model$fullpars[12]/(kap*lam1))*(model$fullpars[13]/(kap*lam2))),
                             model$fullpars[14]*sqrt((model$fullpars[12]/(kap*lam1))*(model$fullpars[13]/(kap*lam2))), model$fullpars[13]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[3]] <- matrix(c(model$fullpars[15]/(kap*lam1), model$fullpars[17]*sqrt((model$fullpars[15]/(kap*lam1))*(model$fullpars[16]/(kap*lam2))),
                             model$fullpars[17]*sqrt((model$fullpars[15]/(kap*lam1))*(model$fullpars[16]/(kap*lam2))), model$fullpars[16]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     covmat[[4]] <- matrix(c(model$fullpars[18]/(kap*lam1), model$fullpars[20]*sqrt((model$fullpars[18]/(kap*lam1))*(model$fullpars[19]/(kap*lam2))),
                             model$fullpars[20]*sqrt((model$fullpars[18]/(kap*lam1))*(model$fullpars[19]/(kap*lam2))), model$fullpars[19]/(kap*lam2)),
                           2,2,byrow=TRUE)
     
     # get decision bound parameters for this subject
     b <- matrix(c(1, model$fullpars[20+(sub-1)*8+3], 
                   1, model$fullpars[20+(sub-1)*8+5],
                   model$fullpars[20+(sub-1)*8+7], 1),
                  3,2,byrow=TRUE)
     c <- matrix(c(model$fullpars[20+(sub-1)*8+4],
                   model$fullpars[20+(sub-1)*8+6],
                   model$fullpars[20+(sub-1)*8+8]),3,1)
     
     if (out_format=="vector"){
       tryCatch(
         {
         predicted <- c(predicted, as.vector(matrix_predict_r3x2(means,covmat,b,c)))
         },
         error=function(e) {
           warning("The fitted GRT-wIND model could not produce predictions")
           predicted <- NA
           }
         )
       } else if (out_format=="list"){
       predicted <- c(predicted, list(matrix_predict_r3x2(means,covmat,b,c)))
     }
   }

   return(predicted)
   
}
