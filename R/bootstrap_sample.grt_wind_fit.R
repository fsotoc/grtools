
#' Obtain bootstrap samples for parameter estimates of a fitted GRT-wIND model
#'
#' Produces N bootstrap samples for parameter estimates of a previously-fitted GRT-wIND model,
#' and returns all samples of parameters and (1-alpha)*100% bootstrap confidence intervals for them.
#' 
#' @param x An object returned by \code{grt_wind_fit} or
#'   \code{grt_wind_fit_parallel}
#' @param N_rows Number of trials per row; it can be a single number (equal number of
#' trials per row; same for all participants), a vector with four numbers (unequal number 
#' of trials per row, but same across participants), or a list (one element per participant)
#' of such single values or four-element vectors. If you have an original dataset in the form
#' of a \code{cmats} object (list of confusion matrices, see \code{grt_wind_fit}), you
#' can compute row counts using \code{lapply(cmats, rowSums)} and then enter the result here.
#' @param N Number of bootstrap samples to be obtained. Defaults to 1,000.
#' @param alpha Error probability parameter used to create (1-alpha)*100\% two-sided confidence
#' interval. Defaults to 0.05 (95\% confidence interval).
#' @param n_cores Number of cores that the program would use to run the bootstrapping, by default
#' it uses all available cores minus one.
#' @param ... Additional arguments passed to methods.
#' 
#' @return An object of class "\code{model_samples}" including all samples and confidence intervals.
#' 
#' @details This procedure uses parametric bootstrap sampling from the fitted GRT-wIND model
#' to build confidence intervals. Each sample of parameters is obtained by sampling random data
#' from the fitted model, fitting a GRT-wIND model to that generated data using \code{grt_wind_fit},
#' and storing the resulting parameter estimates. Confidence intervals are obtained through the function
#' \code{confidence_interval}, which uses a simple percentile method.
#' 
#' @references
#' 
#' @examples 
#' # Create list with confusion matrices. In this example, we enter data from
#' # an experiment with 5 participants. For each participant, inside the c(...),
#' # enter the data from row 1 in the matrix, then from row 2, etc.
#' cmats <- list(matrix(c(161,41,66,32,24,147,64,65,37,41,179,43,14,62,22,202),nrow=4,ncol=4,byrow=TRUE))
#' cmats[[2]] <- matrix(c(126,82,67,25,8,188,54,50,34,75,172,19,7,103,14,176),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[3]] <- matrix(c(117,64,89,30,11,186,69,34,21,81,176,22,10,98,30,162),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[4]] <- matrix(c(168,57,47,28,15,203,33,49,58,54,156,32,9,96,9,186),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[5]] <- matrix(c(169,53,53,25,34,168,69,29,38,48,180,34,19,44,60,177),nrow=4,ncol=4,byrow=TRUE)
#' 
#' # fit the GRT_wIND model to data
#' fitted_model <- grt_wind_fit(cmats)
#' 
#' # get number of rows in original matrices
#' N_rows <- lapply(cmats, rowSums)
#' 
#' # produce 10 bootstrap samples
#' bootstrap_samples <- bootstrap_sample(fitted_model, N_rows=N_rows, N=10, n_cores=1)
#' 
#' # look at the confidence intervals
#' bootstrap_samples$cis
#' 
#' @export
bootstrap_sample <- function(x, ...) UseMethod("bootstrap_sample")


#' @export
bootstrap_sample.grt_wind_fit <- function(x, N_rows, N = 1000, alpha = 0.05, n_cores = 0, ...){
  fitted_model <- x
  # code by Sanaz Hosseini
  if (length(N) != 1 || is.na(N) || N < 1) {
    stop("'N' must be a positive integer")
  }
  N <- as.integer(N)
  
  # Determine how many cores we are going to use
  available_cores <- parallel::detectCores()  #counting available cores of the system
  if (is.na(available_cores) || available_cores < 2) {
    available_cores <- 2
  }
  if (n_cores == 0 || n_cores > available_cores) {  #if no input or input is larger than available cores
    n_cores <- available_cores - 1                  #use all the available cores minus 1 (just to be able to do sth else with the system, Googling stuff :D)
  }
  if (is.na(n_cores) || n_cores < 1) {
    n_cores <- 1
  }
  
  # Run serially if requested to avoid parallel backend issues
  if (n_cores == 1) {
    result <- bsample(fitted_model = fitted_model, N_rows = N_rows, N = N)
    result_ci <- confidence_interval(result, mle_pars = fitted_model$fullpars, alpha = alpha)
    final_output <- list(samples = result, cis = result_ci)
    class(final_output) <- "model_samples"
    return(final_output)
  }
  
  R<-0   # FLAG: for when N is not divisible by n_cores
  
  if (N < n_cores){  # if the number of calculations is smaller than the cores,
    n_cores <- N     # just assign as many cores to calculations as the number of required calculation
    N1<-1            # N1 is the number of samples obtained per core through the function bsample. 
                     # Here 1 because each core will get only one sample.
    
  }
  else{                   # if the number of calculations is not smaller than the cores, 
    if (N%%n_cores==0){   # if the number of calculations is divisible by the number of cores,
      N1<-N/n_cores       # assign each core "N1=number_of_computations/n_cores" tasks to run 
    }
    else{                  # if the number of calculations is not divisible by the number of cores,
      R<-N%%n_cores        # Change the R (FLAG) to the residual of the devision (so that we won't miss these calculations)
      N1<-floor(N/n_cores) # assign each core "N1=number_of_computations/n_cores" tasks to run   
    }
  }
  
  cl <- parallel::makeCluster(n_cores)  # making a cluster with n_cores nodes
  parallel::clusterEvalQ(cl, library(grtools))
  
  PBS <- parallel::clusterCall(cl, fun=bsample, fitted_model=fitted_model , 
                               N=N1, N_rows=N_rows)  # sample using function bsample
  
  if (R!=0){                         # If the FLAG is on, 
    
    cl2 <- parallel::makeCluster(R)   # make a cluster, with smaller number of nodes
    parallel::clusterEvalQ(cl2, library(grtools))
    
    PBS1 <- parallel::clusterCall(cl2, fun=bsample, fitted_model=fitted_model , 
                                  N=1, N_rows=N_rows) #sample using function bsample
    parallel::stopCluster(cl2)
  }
  
  parallel::stopCluster(cl)   #Stop the cluster after finishing
  
  #Producing the output matrix
  
  if (R==0){     #when there is no FLAG
    # create N (number of samples) x K (number of parameters) matrix and populate it with sampled values
    result<-matrix(nrow = N, ncol=ncol(PBS[[1]]), byrow = TRUE)  
    for (i in seq(0, length(PBS)*N1-1, N1)){    
      for (k in 1:N1){                          
        result[i+k, ] <- PBS[[i/N1+1]][k,]
      }
    } 
  }
  else{         #when the flag is on
    # create N (number of samples) x K (number of parameters) matrix and populate it with sampled values
    result<-matrix(nrow = N, ncol=ncol(PBS[[1]]), byrow = TRUE)
    for (i in seq(0, length(PBS)*N1-1, N1)){
      for (k in 1:N1){
        result[i+k, ] <- PBS[[i/N1+1]][k,]
      }
    } 
    for (i in 1:R){
      result[length(PBS)*N1+i, ] <- PBS1[[i]]
    }
  }
  
  # names of the parameters for the group model
  paramether1 <- c("A1B1_mean_A", "A1B1_mean_B", "A2B1_mean_A", "A2B1_mean_B", 
                   "A1B2_mean_A", "A1B2_mean_B", "A2B2_mean_A", "A2B2_mean_B",
                   "A1B1_var_A", "A1B1_var_B", "A1B1_corr","A2B1_var_A", "A2B1_var_B", "A2B1_corr",
                   "A1B2_var_A", "A1B2_var_B", "A1B2_corr","A2B2_var_A", "A2B2_var_B", "A2B2_corr")
  
  # name individual model parameters
  paramether2 <- matrix(vector (mode="character",length=NCOL(result)-20), ncol = 1, byrow = TRUE)
  J<-c()
  # parameters defined for each subject
  J[1]="_kappa"
  J[2]="_lambda"
  J[3]="_by1"
  J[4]="_a1"
  J[5]="_bx2"
  J[6]="_a2"
  J[7]="_by2"
  J[8]="_bx3"
  J[9]="_a3"
  # number of individual parameters is different for r3x2 model (type-2 GRTapas)
  if (any(fitted_model$restrictions=="r3x2")){
    paramether2 <- paste0("subject",floor((row(paramether2)-1)/8)+1, J[rep(c(1:4,7,6,8,9), fitted_model$N)])
  } else {
    paramether2 <- paste0("subject",floor((row(paramether2)-1)/6)+1, J[(row(paramether2)-1)%%6+1])
  }
  paramether <- append(paramether1, paramether2)
  colnames(result) <- paramether
  
  result_ci <- confidence_interval(result, mle_pars=fitted_model$fullpars, alpha=alpha)
  
  final_output <- list(samples=result, cis=result_ci) # return the result matrix as output 
  
  # set class of the output as "model_samples"
  class(final_output) <- "model_samples"
  
  return(final_output)
}
