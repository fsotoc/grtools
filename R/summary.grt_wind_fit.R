#' @export
summary.grt_wind_fit <- function(fitted_model) {

  if (fitted_model$convergence==0){
    cat("The optimization algorithm was successful\n\n")
  } else {
    cat("The optimization algorithm may have failed\n")
    cat("The following message was produced by the optim() function:\n")
    cat(paste("\t",fitted_model$message,"\n\n"))
  }
  
  cat("Measures of Fit: \n")
  cat(paste("\tLog-likelihood:", round(-fitted_model$value, 2), "\n"))
  cat(paste("\tR-squared:", round(fitted_model$R2, 4), "\n"))
  
  if (!is.null(fitted_model$lr_test)){
    cat("\nResults of the Likelihood Ratio tests:\n")
    forprint <- fitted_model$lr_test
    forprint$pval <- round(forprint$pval, digits = 3)
    print(forprint, digits=3, row.names=F)
  }
  
  if (!is.null(fitted_model$wald)){
    cat("\nResults of the Wald tests:\n")
    print(fitted_model$wald_test, digits=3, row.names=F)
  }
  
}