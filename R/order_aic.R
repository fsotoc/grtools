
#'  Compute AIC for a list of models and rank them according to that measure
#' 
#' @param fitted_models A list of fitted GRT or GRT-wIND models objects.
#' @param model_names An array of names for the fitted models
#' @param n_data Total number of unique data points used to fit the models. This
#' is 16 for traditional 2x2 GRT models, 16 times the number of participants
#' for the 2x2 GRT-wIND model, etc.
#' 
#' @export
#' 
order_aic <-function(fitted_models, model_names, n_data) {
  # number of models
  n_mods <- length(fitted_models)
  
  aic_list <- rep(0,n_mods)
  L <- rep(0,n_mods)
  for(i in 1:n_mods){
    L[i] <- -fitted_models[[i]]$value
    m <- length(fitted_models[[i]]$par)
    aic_list[i] <- -2*L[i]+2*m+(2*m^2+2*m)/(n_data-m-1)
  }
  
  aic_exp <- rep(0,n_mods)    
  aic_weight <- rep(0,n_mods)
  aic_exp <- exp(-(aic_list-min(aic_list))/2)
  aic_weight <- aic_exp/sum(aic_exp)
  
  ordered_aic <- data.frame(model_names,L,aic_list,aic_weight)
  colnames(ordered_aic) <- c("model","log-likelihood", "AIC", "AIC weight")
  ordered_aic <- ordered_aic[order(-aic_weight),]
  ordered_aic[4] <- prettyNum(round(ordered_aic[4], digits=3)[[1]], nsmall=2)
  return(ordered_aic)
}
