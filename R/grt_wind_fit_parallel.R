#' Fit a GRT-wIND model to data several times in parallel
#' 
#' It fits a GRT-wIND model to data by running \code{\link{grt_wind_fit}} repeated times,
#' each time with a different value for the starting parameters. It returns the model
#' with a highest log-likelihood from all the runs.
#' 
#' @param n_reps Number of times the optimization algorithm should be run. Must be provided.
#' @param n_cores Number of cores to be used. It defaults to all available cores
#' minus one.
#' @inheritParams grt_wind_fit
#' @details \code{rand_pert} must be higher than zero. For more details, see 
#' \code{\link{grt_wind_fit}}
#' 
#' @return An object of class "\code{grt_wind_fit}."
#' For more details and examples, see \code{\link{grt_wind_fit}}
#' 
#' @references Soto, F. A., Musgrave, R., Vucovich, L., & Ashby, F. G. (2015).
#' General recognition theory with individual differences: A new method for
#' examining perceptual and decisional interactions with an application to
#' face perception. \emph{Psychonomic Bulletin & Review, 22}(1), 88-111.
#' 
#' @export
grt_wind_fit_parallel <- function(cmats, start_params=c(), model="full", rand_pert=0.3,
                                  control=list(),
                                  n_reps, n_cores=0) {
  
  
  # Calculate the number of available cores
  available_cores <- parallel::detectCores()
  
  # if the number of cores was not provided or the user asked for more
  # cores than the available number, then use available minus one
  if ( n_cores == 0 || n_cores > available_cores ) {
    n_cores <- available_cores - 1
  }
  
  # if the number of repetitions is less than the number of cores, set the latter to the former
  if ( n_reps < n_cores ) {
    n_cores <- n_reps
  }
  
  # Initiate cluster and load package in the workers
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, library(grtools))
  
  # Create a list with repetitions of the data
  rep_cmats <- rep(list(cmats), n_reps)
  
  # run in parallel
  results <- parallel::parLapply(cl, rep_cmats, grt_wind_fit, 
                                 start_params=start_params,
                                 model = model,
                                 rand_pert=rand_pert,
                                 control=control)
  
  # stop the cluster
  parallel::stopCluster(cl)
  
  # get the best-fitting model
  nlls <- c()
  for (i in 1:n_reps){
    nlls <- c(nlls, results[[i]]$value)
  }
  best_model <- results[[which.min(nlls)]]
  
  return(best_model)
  
}