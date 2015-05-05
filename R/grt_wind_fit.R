#' Fit a GRT-wIND model to identification data
#' 
#' Uses the BFGS optimization method to fit a full GRT-wIND model to data from a
#' 2x2 identification experiment (see Soto et al., 2015).
#' 
#' @param cmats List of confusion matrices. Each entry in the list should contain
#'   the 4x4 confusion matrix from one individual (see Details).
#' @param start_params An optional vector of parameters to start the
#'   optimization algorithm. You can provide either the group parameters or both
#'   group and individual parameters. It will be created automatically if not
#'   provided (see Details).
#' @param rand_pert Maximum value of a random perturbation added to the starting
#'   parameters. With a value of zero, the algorithm is started exactly at
#'   \code{start_params}. As the value of \code{rand_pert} is increased, the starting
#'   parameters become closer to be "truly random."
#' @param control A list of control parameters for the \code{optim} function. See
#'   \code{\link[stats]{optim}}.
#' 
#' @return An object of class "\code{grt_wind_fit}."
#' 
#' The function \code{\link{wald}} is used to obtain the
#' results of Wald tests of perceptual separability, decisional separability and
#' perceptual independence, using the maximum likelihood parameters estimated
#' through \code{grt_wind_fit}.
#' 
#' The function \code{summary} is used to obtain a summary of the model fit to 
#' data and the results of Wald tests (if available). The function
#' \code{\link[=plot.grt_wind_fit]{plot}} is used to print a graphical
#' representation of the fitted model.
#' 
#' @details 
#' We recommend that you fit GRT-wIND to your data repeated times, each time 
#' with a different value for the starting parameters, and keep the solution
#' with the smallest negative log-likelihood. This facilitates finding the true 
#' maximum likelihood estimates of the model's parameters, which is necessary to
#' make valid conclusions about dimensional separability and independence. The 
#' function \code{\link{grt_wind_fit_parallel}} does exactly this, taking 
#' advantage of multiple CPUs in your computer to speed up the process. For most
#' applications, you should use \code{\link{grt_wind_fit_parallel}} instead of
#' the current function.
#' 
#' A 2x2 identification experiment involves two dimensions, A and B, each with
#' two levels, 1 and 2. Stimuli are represented by their level in each dimension
#' (A1B1, A1B2, A2B1, and A2B2) and so are their corresponding correct
#' identification responses (a1b1, a1b2, a2b1, and a2b2).
#' 
#' The data from a single participant in the experiment should be ordered in a
#' 4x4 confusion matrix with rows representing stimuli and columns representing
#' responses. Each cell has the frequency of responses for the stimulus/response
#' pair. Rows and columns should be ordered in the following way:
#' 
#' \itemize{ \item{Row 1: Stimulus A1B1} \item{Row 2: Stimulus A2B1} 
#' \item{Row 3: Stimulus A1B2} \item{Row 4: Stimulus A2B2} \item{Column
#' 1: Response a1b1} \item{Column 2: Response a2b1} \item{Column 3: Response a1b2} 
#' \item{Column 4: Response a2b2} }
#' 
#' The argument \code{cmats} is a list with all individual confusion matrices
#' from an experimental group.
#' 
#' If the value of \code{start_params} is not provided, the starting parameters
#' for the optimization algorithm are the following: 
#' \itemize{ 
#' \item{Means:}{ A1B1=(0,0), A2B1=(1,0), A1B2=(1,0), A2B1=(1,1)} 
#' \item{Variances:}{ All set to one} 
#' \item{Correlations:}{ All set to zero} 
#' \item{Decision bounds:}{ All orthogonal to their corresponding dimension} 
#' \item{Attention parameters:}{ No differences in global attention (kappa=2) and no 
#' selective attention (lambda=0.5)}
#' }
#' 
#' Note that a random value will be added to these parameters if
#' \code{rand_pert} is given a value higher than zero.
#' 
#' @references Soto, F. A., Musgrave, R., Vucovich, L., & Ashby, F. G. (2015).
#'   General recognition theory with individual differences: A new method for
#'   examining perceptual and decisional interactions with an application to
#'   face perception. \emph{Psychonomic Bulletin & Review, 22}(1), 88-111.
#'    
#' @examples 
#' # Create list with confusion matrices # In this example, we enter data from
#' # an experiment with 5 participants. For each participant, inside the c(...),
#' # enter the data from row 1 in the matrix, then from row 2, etc.
#' cmats <- list(matrix(c(100,1,9,8,10,110,7,4,31,3,80,10,54,4,52,19),nrow=4,ncol=4,byrow=TRUE))
#' cmats[[2]] <- matrix(c(122,7,0,1,1,102,1,5,3,2,111,9,11,7,11,106),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[3]] <- matrix(c(107,0,5,0,0,101,1,4,3,1,113,2,0,3,1,108),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[4]] <- matrix(c(122,1,0,0,1,120,0,0,0,0,118,6,0,1,6,118),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[5]] <- matrix(c(89,17,6,4,4,81,8,6,14,7,86,1,11,25,17,26),nrow=4,ncol=4,byrow=TRUE)
#' 
#' # fit the model to data
#' fitted_model <- grt_wind_fit(cmats)
#' 
#' # plot a graphical representation of the fitted model
#' plot(fitted_model)
#' 
#' # optionally, you can run a Wald test of separability and independence
#' fitted_model <- wald(fitted_model, cmats)
#' 
#' # print a summary of the fitted model and tests to screen
#' summary(fitted_model)
#' 
#' # a data frame with the values of all individual parameters is available in
#' # in fitted_model$indpar
#' fitted_model$indpar
#' 
#' @export
grt_wind_fit <- function(cmats, start_params=c(), rand_pert=0.3, control=list()) {  
  
  # get number of subjects
  N = length(cmats)
  

  #---------------------------------------------------------
  # define upper and lower limits for the parameters
  mu_up <- Inf
  mu_low <- -Inf
  corr_up <- 1
  corr_low <- -1
  var_up <- Inf
  var_low <- 0
  lambda_low <- 0.001
  lambda_up <- 0.999
  kappa_low <- 0
  kappa_up = Inf
  b_up <- Inf
  b_low <- -Inf
  a_up <- Inf
  a_low <- -Inf
  
  #----------------------------------------------------------
  # get starting values if they were not provided in start_params
  
  if (length(start_params)==0){
    
    # choose values for group parameters
    start_params = c(1,0,0,1,1,1,     # all means either 0 or 1 (PS)
                     0,               # A1B1 correlation (PI)
                     1,1,             # A1B2 variances (PS)
                     0,               # A1B2 correlation (PI)
                     1,1,             # A2B1 variances (PS)
                     0,               # A2B1 correlation (PI)
                     1,1,             # A2B2 variances (PS)
                     0)               # A2B2 correlation (PI)
    
    # choose values for individual parameters
    lambda_start <- rep(0.5,N)                                 # start with no selective attention
    kappa_start <- rep(2,N)                                    # start with no global attention
    bx_start <- rep(0,N)                                # start with bounds orthogonal to dimension (DS)
    by_start <- rep(0,N)
    ax_start <- rep(-0.5,N)                             # start with bound right between the two means
    ay_start <- rep(-0.5,N)
    
  # if only guesses for the group parameters were provided,
  # then initialize the individual parameters
  } else if (length(start_params)==16) {
    
    lambda_start <- rep(0.5,N)                        # start with no selective attention
    kappa_start <- rep(2,N)                           # start with no global attention
    bx_start <- rep(0,N)                              # start with bounds orthogonal to dimension (DS)
    by_start <- rep(0,N)
    ax_start <- rep(-0.5,N)                           # start with bound right between the two means
    ay_start <- rep(-0.5,N)
    
  }
  
  # lower and upper limits for group parameters
  low_params = c(mu_low, mu_low, mu_low, mu_low, mu_low, mu_low, corr_low, var_low, var_low, corr_low, var_low, var_low, corr_low, var_low, var_low, corr_low)
  up_params = c(mu_up, mu_up, mu_up, mu_up, mu_up, mu_up, corr_up, var_up, var_up, corr_up, var_up, var_up, corr_up, var_up, var_up, corr_up)
  
  
  # now that all possibilities are covered, put together vector with all
  # starting parameters and limits
  for (i in 1:N) {
    start_params <- c(start_params, kappa_start[i], lambda_start[i], bx_start[i], ax_start[i], by_start[i], ay_start[i])
    low_params <- c(low_params, kappa_low, lambda_low, b_low, a_low, b_low, a_low)
    up_params <- c(up_params, kappa_up, lambda_up, b_up, a_up, b_up, a_up)
  }

  # add random perturbations to parameters with maximum value rand_pert
  start_params <- start_params+(runif(length(start_params))*(2*rand_pert)-rand_pert)
  
  
  #--------------------------------------------------------
  # set up default control parameters for optimization
  # these values have yielded the highest likelihood in problems
  # using real data, but this might be specific to each data set
  if ( !("maxit" %in% names(control)) ){
    control$maxit <- 1e+5
  }
  if ( !("ndeps" %in% names(control)) ){
    control$ndeps <- rep(1e-2, times=length(start_params))
  }
  if ( !("factr" %in% names(control)) ){
    control$factr <- 1e-5
  }
  
  # find maximum likelihood estimates
  results <- optim(start_params, grt_wind_nll, data=cmats, method='L-BFGS-B', lower=low_params, upper=up_params, control=control, hessian = F)
  
  # put in a nice list
  results$means <- matrix(c(0,0,results$par[1:6]), nrow=4, ncol=2, byrow=T)
  results$covmat <- list()
  results$covmat[[1]] <- matrix(c(1, results$par[7], results$par[7], 1),2,2,byrow=TRUE)
  results$covmat[[2]] <- matrix(c(results$par[8], results$par[10]*sqrt(results$par[8]*results$par[9]),
                          results$par[10]*sqrt((results$par[8])*(results$par[9])), results$par[9]),
                        2,2,byrow=TRUE)
  results$covmat[[3]] <- matrix(c(results$par[11], results$par[13]*sqrt(results$par[11]*results$par[12]),
                          results$par[13]*sqrt(results$par[11]*results$par[12]), results$par[12]),
                        2,2,byrow=TRUE)
  results$covmat[[4]] <- matrix(c(results$par[14], results$par[16]*sqrt(results$par[14]*results$par[15]),
                          results$par[16]*sqrt(results$par[14]*results$par[15]), results$par[15]),
                        2,2,byrow=TRUE)
  results$corr <- results$par[c(7,10,13,16)]
  results$var <- matrix(c(1,1,results$par[c(8,9,11,12,14,15)]), nrow=4, ncol=2, byrow=T)
  
  # individual parameters
  results$indpar <- data.frame()
  pcount <- 16
  for (i in 1:N){
    results$indpar <- rbind(results$indpar,
                            data.frame(subject=i,
                              kappa=results$par[pcount+1],
                              lambda=results$par[pcount+2],
                              bx1=1,
                              by1=results$par[pcount+3],
                              a1=results$par[pcount+4],
                              bx2=results$par[pcount+5],
                              by2=1,
                              a2=results$par[pcount+6]) )
    pcount <- pcount+6
  }
  
  # make the results object class "grt_wind_fit"
  class(results) <- "grt_wind_fit"
  
  results$N <- N
  results$model <- "GRT-wIND"
  results$predicted <- fitted(results)
  results$observed <- c()
  
  for (i in 1:N){
    results$observed <- c(results$observed, as.vector(pmatrix(cmats[[i]])))    
  }
  
  results$R2 <- cor(results$predicted, results$observed)^2
  
  
  return(results)
  
}




