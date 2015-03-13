
#' Fit a hierarchy of traditional GRT models to data
#' 
#' Fits a hierarchy of traditional GRT models to data from a 2x2 identification 
#' experiment, using the BFGS optimization method (See Ashby & Soto, 2015). It
#' then selects the best-fitting model using the AIC. 
#' 
#' @param data A 4x4 confusion matrix (see Details).
#' @param rand_pert Maximum value of a random perturbation added to the starting
#'   parameters. Defaults to 0.3. With a value of zero, the optimization is started exactly at the 
#'   default starting parameters (see Details). As the value of rand_pert is 
#'   increased, the starting parameters become closer to be "truly random."
#' @param n_reps Number of times the optimization algorithm should be run, each time
#' with a different value for the starting parameters. The function will return the
#' model with a highest log-likelihood from all the runs. The value of \code{nreps}
#' defaults to five.
#' @param control A list of optional control parameters for the \code{optim} function. See 
#'   \code{\link[stats]{optim}}. Note that the parameter \code{ndeps} entered 
#'   here should be a single number instead of the vector that is usually passed 
#'   to \code{optim}. This single value is repeated inside \code{grt_hm_fit} to 
#'   create the appropriate vectors.
#'   
#' @return An object of class "\code{grt_hm_fit}."
#'   
#'   The function \code{summary} is used to obtain a summary of results from the
#'   model fit and selection process, including the best-fitting model and
#'   conclusions about perceptual separability and perceptual independence
#'   (decisional separability is assumed by all models)
#'   
#'   The function \code{\link[=plot.grt_hm_fit]{plot}} is used to print a
#'   graphical representation of the best-fitting model.
#'   
#' @details A 2x2 identification experiment involves two dimensions, A and B,
#' each with two levels, 1 and 2. Stimuli are represented by their level in each
#' dimension (A1B1, A1B2, A2B1, and A2B2) and so are their corresponding correct
#' identification responses (a1b1, a1b2, a2b1, and a2b2).
#' 
#' The data from a single participant in the experiment should be ordered in a 
#' 4x4 confusion matrix with rows representing stimuli and columns representing 
#' responses. Each cell has the frequency of responses for the stimulus/response
#' pair. Rows and columns should be ordered in the following way:
#' 
#' \itemize{ \item{Stimulus/Row 1: A1B1} \item{Stimulus/Row 2: A2B1} 
#' \item{Stimulus/Row 3: A1B2} \item{Stimulus/Row 4: A2B2} \item{Response/Column
#' 1: a1b1} \item{Response/Column 2: a2b1} \item{Response/Column 3: a1b2} 
#' \item{Respones/Column 4: a2b2} }
#' 
#' The default starting parameters for the optimization algorithm are the
#' following: \itemize{ \item{Means:}{ A1B1=(0,0), A2B1=(1,0), A1B2=(1,0),
#' A2B1=(1,1)} \item{Variances:}{ All set to one} \item{Correlations:}{ All set
#' to zero} }
#' 
#' Decisional separability is assumed for all models (i.e., decision bounds are
#' fixed and orthogonal to the dimension they divide)
#' 
#' Note that a random value will be added to the default starting parameters if 
#' \code{rand_pert} is given a value higher than zero.
#' 
#' @references Ashby, F. G., & Soto, F. A. (2015). Multidimensional signal
#'   detection theory. In J. R. Busemeyer, J. T. Townsend, Z. J. Wang, & A.
#'   Eidels (Eds.), \emph{Oxford handbook of computational and mathematical
#'   psychology} (pp. 13-34). Oxford University Press: New York, NY.
#'   
#' @examples 
#' # Create a confusion matrix
#' # Inside the c(...) below, we enter the data from row 1 in the 
#' matrix, then from row 2, etc.
#' conf_matrix <- matrix(c(100,1,9,8,
#'                        10,110,7,4,
#'                        31,3,80,10,
#'                        54,4,52,19),
#'                        nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Perform model fit and selection
#' hm_fit_results <- grt_hm_fit(conf_matrix)
#' 
#' # See a summary of the fitting and selection results
#' summary(hm_fit_results)
#' 
#' # plot a graphical representation of the best-fitting model
#' plot(hm_fit_results)
#' 
#' @export
#' 
grt_hm_fit <- function(data, rand_pert=0.3, nreps=5, control=list()){
  
  # fit all models 
  fitted_models <- fit_grt_models(data, rand_pert=rand_pert, nreps=nreps, control=control)
  
  # order all models according to AIC
  o_aic <- order_aic(fitted_models)
  
  # put together a list with best-model parameters for output
  best_model <- list()
  best_model$means <- matrix(0, 4, 2, byrow=TRUE)
  best_model$covmat <- list()
  best_model$a1 <- 0
  best_model$a2 <- 0
  best_model$model <- ""
  
  o_match <- o_aic[1,]$model
  best_model$model <- paste("GRT-", o_match, sep="")
  model_list <- c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", 
                  "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", "{PI, DS}", 
                  "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", 
                  "{1_RHO, DS}", "{PS(B), DS}", "{DS}")
  model_num <- pmatch(o_match, model_list)
  w <- fitted_models[[model_num]]$par
  best_model$model <- paste("GRT-", o_match, sep="")
  switch(model_num,
         mod1={best_model$means[2,1] <- w[1]
               best_model$means[3,2] <- w[2]
               best_model$means[4,1] <- w[1]
               best_model$means[4,2] <- w[2]
               best_model$covmat[[1]] <- diag(2)
               best_model$covmat[[2]] <- diag(2)
               best_model$covmat[[3]] <- diag(2)
               best_model$covmat[[4]] <- diag(2)
               best_model$a1 <- w[3]
               best_model$a2 <- w[4]},
         
         mod2={best_model$means[2,1] <- w[1]
               best_model$means[2,2] <- w[2]
               best_model$means[3,2] <- w[3]
               best_model$means[4,1] <- w[1]
               best_model$means[4,2] <- w[4]
               best_model$covmat[[1]] <- diag(2)
               best_model$covmat[[2]] <- diag(2)
               best_model$covmat[[3]] <- diag(2)
               best_model$covmat[[4]] <- diag(2)
               best_model$a1 <- w[5]
               best_model$a2 <- w[6]},
         
         mod3={best_model$means[2,1] <- w[1]
               best_model$means[3,1] <- w[2]
               best_model$means[3,2] <- w[3]
               best_model$means[4,1] <- w[4]
               best_model$means[4,2] <- w[3]
               best_model$covmat[[1]] <- diag(2)
               best_model$covmat[[2]] <- diag(2)
               best_model$covmat[[3]] <- diag(2)
               best_model$covmat[[4]] <- diag(2)
               best_model$a1 <- w[5]
               best_model$a2 <- w[6]},
         
         mod4={best_model$means[2,1] <- w[1]
               best_model$means[3,2] <- w[2]
               best_model$means[4,1] <- w[1]
               best_model$means[4,2] <- w[2]
               vals <- c(1,w[3],w[3],1)
               best_model$covmat[[1]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[2]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[3]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[4]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$a1 <- w[4]
               best_model$a2 <- w[5]},
         
         
         mod5={best_model$means[2,1] <- w[1]
               best_model$means[2,2] <- w[2]
               best_model$means[3,2] <- w[3]
               best_model$means[4,1] <- w[1]
               best_model$means[4,2] <- w[4]
               vals <- c(1,w[5],w[5],1)
               best_model$covmat[[1]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[2]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[3]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[4]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$a1 <- w[6]
               best_model$a2 <- w[7]},
         
         mod6={best_model$means[2,1] <- w[1]
               best_model$means[2,2] <- w[2]
               best_model$means[3,1] <- w[3]
               best_model$means[3,2] <- w[4]
               best_model$means[4,1] <- w[5]
               best_model$means[4,2] <- w[6]
               best_model$covmat[[1]] <- diag(2)
               best_model$covmat[[2]] <- diag(2)
               best_model$covmat[[3]] <- diag(2)
               best_model$covmat[[4]] <- diag(2)
               best_model$a1 <- w[7]
               best_model$a2 <- w[8]},
         
         mod7={best_model$means[2,1] <- w[1]
               best_model$means[3,1] <- w[2]
               best_model$means[3,2] <- w[3]
               best_model$means[4,1] <- w[4]
               best_model$means[4,2] <- w[3]
               vals <- c(1,w[5],w[5],1)
               best_model$covmat[[1]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[2]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[3]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$covmat[[4]] <- matrix(vals, 2, 2, byrow=TRUE)
               best_model$a1 <- w[6]
               best_model$a2 <- w[7]},
         
         mod8={best_model$means[2,1] <- w[1]
               best_model$means[3,2] <- w[2]
               best_model$means[4,1] <- w[1]
               best_model$means[4,2] <- w[2]
               best_model$covmat[[1]] <- matrix(c(1, w[3], w[3], 1), 2, 2, byrow=TRUE)
               best_model$covmat[[2]] <- matrix(c(1, w[4], w[4], 1), 2, 2, byrow=TRUE)
               best_model$covmat[[3]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow=TRUE)
               best_model$covmat[[4]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow=TRUE)
               best_model$a1 <- w[7]
               best_model$a2 <- w[8]},
         
         mod9={best_model$means[2,1] <- w[1]
               best_model$means[2,2] <- w[2]
               best_model$means[3,2] <- w[3]
               best_model$means[4,1] <- w[1]
               best_model$means[4,2] <- w[4]
               best_model$covmat[[1]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow=TRUE)
               best_model$covmat[[2]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow=TRUE)
               best_model$covmat[[3]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow=TRUE)
               best_model$covmat[[4]] <- matrix(c(1, w[8], w[8], 1), 2, 2, byrow=TRUE)
               best_model$a1 <- w[9]
               best_model$a2 <- w[10]},
         
         mod10={best_model$means[2,1] <- w[1]
                best_model$means[2,2] <- w[2]
                best_model$means[3,1] <- w[3]
                best_model$means[3,2] <- w[4]
                best_model$means[4,1] <- w[5]
                best_model$means[4,2] <- w[6]
                vals <- c(1,w[7],w[7],1)
                best_model$covmat[[1]] <- matrix(vals, 2, 2, byrow=TRUE)
                best_model$covmat[[2]] <- matrix(vals, 2, 2, byrow=TRUE)
                best_model$covmat[[3]] <- matrix(vals, 2, 2, byrow=TRUE)
                best_model$covmat[[4]] <- matrix(vals, 2, 2, byrow=TRUE)
                best_model$a1 <- w[8]
                best_model$a2 <- w[9]},
         
         mod11={best_model$means[2,1] <- w[1]
                best_model$means[3,1] <- w[2]
                best_model$means[3,2] <- w[3]
                best_model$means[4,1] <- w[4]
                best_model$means[4,2] <- w[3]
                best_model$covmat[[1]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow=TRUE)
                best_model$covmat[[2]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow=TRUE)
                best_model$covmat[[3]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow=TRUE)
                best_model$covmat[[4]] <- matrix(c(1, w[8], w[8], 1), 2, 2, byrow=TRUE)
                best_model$a1 <- w[9]
                best_model$a2 <- w[10]},
         
         mod12={best_model$means[2,1] <- w[1]
                best_model$means[2,2] <- w[2]
                best_model$means[3,1] <- w[3]
                best_model$means[3,2] <- w[4]
                best_model$means[4,1] <- w[5]
                best_model$means[4,2] <- w[6]
                best_model$covmat[[1]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow=TRUE)
                best_model$covmat[[2]] <- matrix(c(1, w[8], w[8], 1), 2, 2, byrow=TRUE)
                best_model$covmat[[3]] <- matrix(c(1, w[9], w[9], 1), 2, 2, byrow=TRUE)
                best_model$covmat[[4]] <- matrix(c(1, w[10], w[10], 1), 2, 2, byrow=TRUE)
                best_model$a1 <- w[11]
                best_model$a2 <- w[12]})
  
  # save convergence info for best-fitting model
  best_model$convergence <- fitted_models[[model_num]]$convergence
  best_model$message <- fitted_models[[model_num]]$message
  
  # get observed and predicted values
  best_model$predicted <- as.vector(matrix_predict(
    best_model$means, best_model$covmat, diag(2), 
    matrix(c(best_model$a1, best_model$a2), 2, 1)) )
  best_model$observed <- as.vector(pmatrix(data))
  
  # return object of class grt_hm_fit
  results <- list(table=o_aic, best_model=best_model)
  class(results) <- "grt_hm_fit"
  return (results)
}

########################################################
# Function that actually performs maximum-likelihood estimation
# for all models:

fit_grt_models <- function(data, rand_pert=0, nreps=1, control=control){     
  
  # if ndeps was not selected by the user, assign a default value
  if (is.null(control[["ndeps"]])) {
    control$ndeps <- 1e-1
  }
  
  
  # fit model 1
  start_params <- c(1, 1, -.5, -.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf, -Inf, -Inf, -Inf)
  high_params <- c(Inf, Inf, Inf, Inf)
  min_nll <- Inf
  for(i in 1:nreps) {
    init_par <- rand_start(start_params, low_params, high_params, rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod1, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model1 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 2
  start_params <- c(1,0,1,1,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod2, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model2 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 3
  start_params <- c(1,0,1,1,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod3, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model3 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 4
  start_params <- c(1,1,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod4, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model4 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 5
  start_params <- c(1,0,1,1,0,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod5, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model5 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 6
  start_params <- c(1,0,0,1,1,1,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod6, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model6 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 7
  start_params <- c(1,0,1,1,0,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod7, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model7 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 8
  start_params <- c(1,1,0,0,0,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,1,1,1,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod8, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model8 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 9
  start_params <- c(1,0,1,1,0,0,0,0,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod9, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model9 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 10
  start_params <- c(1,0,0,1,1,1,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod10, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model10 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 11
  start_params <- c(1,0,1,1,0,0,0,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod11, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model11 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 12
  start_params <- c(1,0,0,1,1,1,0,0,0,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf) 
  min_nll <- Inf
  for(i in 1:nreps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=negloglik_mod12, data=data, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model12 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  fitted_models <- list(mle_model1, mle_model2, mle_model3, mle_model4, mle_model5, 
                     mle_model6, mle_model7, mle_model8, mle_model9, mle_model10,
                     mle_model11, mle_model12)
  return(fitted_models)
}


##############################################
# Function that computes AIC and ranks models according to it:

order_aic <-function(fitted_models) {
  aic_list <- rep(0,12)
  L <- rep(0,12)
  model <- c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", "{PI, DS}", "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", "{1_RHO, DS}", "{PS(B), DS}", "{DS}")
  for(i in 1:12){
    L[i] <- -fitted_models[[i]]$value
    m <- length(fitted_models[[i]]$par)
    aic_list[i] <- -2*L[i]+2*m+(2*m^2+2*m)/(16-m-1)
  }
  
  aic_exp <- rep(0,12)    
  aic_weight <- rep(0,12)
  aic_exp <- exp(-(aic_list-min(aic_list))/2)
  aic_weight <- aic_exp/sum(aic_exp)
  
  ordered_aic <- data.frame(model,L,aic_list,aic_weight)
  colnames(ordered_aic) <- c("model","log-likelihood", "AIC", "AIC weight")
  ordered_aic <- ordered_aic[order(-aic_weight),]
  ordered_aic[4] <- prettyNum(round(ordered_aic[4], digits=3), nsmall=2)
  return(ordered_aic)
}
