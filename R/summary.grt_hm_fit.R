#' @export
summary.grt_hm_fit <- function(object, ...) {
  hm_list <- object
  
  if (hm_list$best_model$convergence==0){
    cat("The optimization algorithm was successful\n\n")
  } else {
    cat("The optimization algorithm may have failed\n")
    cat("The following message was produced by the optim() function:\n")
    cat(paste("\t", hm_list$best_model$message, "\n\n"))
  }
  
  # get check table to print conclusions
  check_table <- data.frame( model=c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", 
                            "{PI, DS}", "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", "{1_RHO, DS}", "{PS(B), DS}", "{DS}") )
  check_table$PS_A <- c("yes", "yes", "no", "yes", "yes", "no", "no", "yes", "yes", "no", "no", "no")
  check_table$PS_B <- c("yes", "no", "yes", "yes", "no", "no", "yes", "yes", "no", "no", "yes", "no")
  check_table$PI <- c("yes", "yes", "yes", "no", "no", "yes", "no", "no", "no", "no", "no", "no" )
  
  cat("Summary of measures of fit for all models:\n")
  cat("(Models are ranked according to AIC)\n\n")
  names(hm_list$table) <- c("Model", "  Log-likelihood", "         AIC", "AIC weight")
  print(hm_list$table, row.names=F)
  
  cat("\n")
  best_model_str <- hm_list$table[1,1]
  cat(paste("Best fitting model:", best_model_str, "\n"))
  cat(paste("Perceptual Separability of A?:", check_table[check_table$model==best_model_str, 2], "\n"))
  cat(paste("Perceptual Separability of B?:", check_table[check_table$model==best_model_str, 3], "\n"))
  cat(paste("Perceptual Independence?:", check_table[check_table$model==best_model_str, 4], "\n"))
  
  
}
