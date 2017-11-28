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
#' @param model A string indicating what GRT-wIND model to run. By default is the
#' full model ("full"), but restricted models are also available. "PS(A)" and "PS(B)"
#' are models that assume PS for dimension A and dimension B, respectively. "PI" is
#' a model that assumes perceptual independence.  "DS(A)" and "DS(B)" are models
#' that assume DS for dimension A ("PS(A)") and dimension B ("PS(B)"), respectively.
#' "1-VAR" is a model that assumes that all variances are equal to one.
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
#' cmats <- list(matrix(c(161,41,66,32,24,147,64,65,37,41,179,43,14,62,22,202),nrow=4,ncol=4,byrow=TRUE))
#' cmats[[2]] <- matrix(c(126,82,67,25,8,188,54,50,34,75,172,19,7,103,14,176),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[3]] <- matrix(c(117,64,89,30,11,186,69,34,21,81,176,22,10,98,30,162),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[4]] <- matrix(c(168,57,47,28,15,203,33,49,58,54,156,32,9,96,9,186),nrow=4,ncol=4,byrow=TRUE)
#' cmats[[5]] <- matrix(c(169,53,53,25,34,168,69,29,38,48,180,34,19,44,60,177),nrow=4,ncol=4,byrow=TRUE)
#' 
#' # fit the model to data
#' fitted_model <- grt_wind_fit(cmats)
#' 
#' # plot a graphical representation of the fitted model
#' plot(fitted_model)
#' 
#' # optionally, you can run likelihood ratio tests of separability and independence
#' # (recommended method, but very slow)
#' fitted_model <- lr_test(fitted_model, cmats)
#' 
#' # alternatively, you can run a Wald test of separability and independence 
#' # (less recommended method, due to common problems with numerical estimation, but faster)
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
grt_wind_fit <- function(cmats, start_params=c(), model="full", rand_pert=0.3, control=list()) {  
  
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
    start_params = c(0,0,1,0,0,1,1,1,     # all means either 0 or 1 (first mean vector fixed + PS)
                     1,1,             # A1B1 variances (fixed)
                     0,               # A1B1 correlation (PI)
                     1,1,             # A1B2 variances (PS)
                     0,               # A1B2 correlation (PI)
                     1,1,             # A2B1 variances (PS)
                     0,               # A2B1 correlation (PI)
                     1,1,             # A2B2 variances (PS)
                     0)               # A2B2 correlation (PI)
    
  } 
  
  # if we only have group parameters, initialize individual parameters
  if (length(start_params)==20) {
    
    lambda_start <- rep(0.5,N)                        # start with no selective attention
    kappa_start <- rep(2,N)                           # start with no global attention
    bx_start <- rep(0,N)                              # start with bounds orthogonal to dimension (DS)
    by_start <- rep(0,N)
    ax_start <- rep(-0.5,N)                           # start with bound right between the two means
    ay_start <- rep(-0.5,N)
    
    for (i in 1:N) {
      start_params <- c(start_params, kappa_start[i], lambda_start[i], bx_start[i], ax_start[i], by_start[i], ay_start[i])
    }
    
  }
  
  #----------------------------------------------------------
  # lower and upper limits
  # group parameters
  low_params = c(mu_low, mu_low, mu_low, mu_low, mu_low, mu_low, mu_low, mu_low,  
                 var_low, var_low, corr_low, 
                 var_low, var_low, corr_low, 
                 var_low, var_low, corr_low, 
                 var_low, var_low, corr_low)
  up_params = c(mu_up, mu_up, mu_up, mu_up, mu_up, mu_up, mu_up, mu_up, 
                var_up, var_up, corr_up, 
                var_up, var_up, corr_up, 
                var_up, var_up, corr_up, 
                var_up, var_up, corr_up)
  
  # individual parameters
  for (i in 1:N) {
    low_params <- c(low_params, kappa_low, lambda_low, b_low, a_low, b_low, a_low)
    up_params <- c(up_params, kappa_up, lambda_up, b_up, a_up, b_up, a_up)
  }
  
  
  # now that all possibilities are covered, put together vector with all
  # starting parameters and limits
  
  
  #--------------------------------------------------------
  # create arrays indicating what parameters are fixed, depending on model
  #
  # 'fixed' is a vector of numerical values.
  # There are three options
  # (a) To indicate that the parameter in fn() is optimised over, use 0
  # (b) To indicate that the parameter in fn() is fixed to its starting value, 
  # use its index. Thus, if the Nth parameter in fn() is fixed to the starting
  # value, then fixed[N] <- N.
  # (c) To indicate that the parameter in fn() is fixed to the value of another
  # parameter, use the index of that parameter. Thus, if the Nth parameter in fn()
  # should have the same value as the Mth parameter, then fixed[N] <- M.
  
  fixed <- rep(0, length(start_params))
  # full
  # for all models, including full, we fix first two means and variances
  fixed[1] <- 1
  fixed[2] <- 2
  fixed[9] <- 9
  fixed[10] <- 10
  
  # PS(A)
  if (model=="PS(A)"){
    # fix means along x axis
    fixed[5] <- 1
    fixed[7] <- 3
    # fix variances along x axis
    fixed[15] <- 9
    fixed[18] <- 12
    
  # PS(B)
  } else if (model=="PS(B)"){
    # fix means along y axis
    fixed[4] <- 2
    fixed[8] <- 6
    # fix variances along y axis
    fixed[13] <- 10
    fixed[19] <- 16
    
  # PI
  } else if (model=="PI") {
    # make sure that the starting value of correlations is zero
    start_params[c(11, 14, 17, 20)] <- 0
    # set value of correlations to starting value
    fixed[c(11, 14, 17, 20)] <- c(11, 14, 17, 20)
    
  # DS(A)
  } else if (model=="DS(A)"){
    for (i in 1:N){
      # make sure that the starting bound is orthogonal
      start_params[20+(i-1)*6+3] <- 0
      # set value of the bound slope to starting value
      fixed[20+(i-1)*6+3] <- 20+(i-1)*6+3
    }
    
  # DS(B)
  } else if (model=="DS(B)"){
    for (i in 1:N){
      # make sure that the starting bound is orthogonal
      start_params[20+(i-1)*6+5] <- 0
      # set value of the bound slope to starting value
      fixed[20+(i-1)*6+5] <- 20+(i-1)*6+5
    }
    
  # 1-VAR
  } else if (model=="1-VAR"){
    # make sure that all starting variances are equal to one
    start_params[c(9,10,12,13,15,16,18,19)]
    # set value of variances to starting value
    fixed[c(9,10,12,13,15,16,18,19)] <- c(9,10,12,13,15,16,18,19)
  }
    
    
    
  # add random perturbations to non-fixed parameters with maximum value rand_pert
  indx = 1:length(start_params)
  start_params[!fixed] <- start_params[!fixed]+(runif(length(start_params[!fixed]))*(2*rand_pert)-rand_pert)
  
  
  #--------------------------------------------------------
  # set up default control parameters for optimization
  # these values have yielded the highest likelihood in problems
  # using real data, but this might be specific to each data set
  if ( !("maxit" %in% names(control)) ){
    control$maxit <- 1e+5
  }
  if ( !("ndeps" %in% names(control)) ){
    control$ndeps <- rep(1e-2, times=sum(fixed==0))
  }
  if ( !("factr" %in% names(control)) ){
    control$factr <- 1e+5
  }
  
  
  
  #--------------------------------------------------------
  # find maximum likelihood estimates
  results <- optifix(start_params, fixed, grt_wind_nll, data=cmats, method='L-BFGS-B', lower=low_params, upper=up_params, control=control, hessian = F)
  
  # put in a nice list
  results$means <- matrix(results$fullpars[1:8], nrow=4, ncol=2, byrow=T)
  results$covmat <- list()
  results$covmat[[1]] <- matrix(c(results$fullpars[9], results$fullpars[11]*sqrt(results$fullpars[9]*results$fullpars[10]),
                                  results$fullpars[11]*sqrt((results$fullpars[9])*(results$fullpars[10])), results$fullpars[9]),
                                2,2,byrow=TRUE)
  results$covmat[[2]] <- matrix(c(results$fullpars[12], results$fullpars[14]*sqrt(results$fullpars[12]*results$fullpars[13]),
                          results$fullpars[14]*sqrt((results$fullpars[12])*(results$fullpars[13])), results$fullpars[13]),
                        2,2,byrow=TRUE)
  results$covmat[[3]] <- matrix(c(results$fullpars[15], results$fullpars[17]*sqrt(results$fullpars[15]*results$fullpars[16]),
                          results$fullpars[17]*sqrt(results$fullpars[15]*results$fullpars[16]), results$fullpars[16]),
                        2,2,byrow=TRUE)
  results$covmat[[4]] <- matrix(c(results$fullpars[18], results$fullpars[20]*sqrt(results$fullpars[18]*results$fullpars[19]),
                          results$fullpars[20]*sqrt(results$fullpars[18]*results$fullpars[19]), results$fullpars[19]),
                        2,2,byrow=TRUE)
  results$corr <- results$fullpars[c(11,14,17,20)]
  results$var <- matrix(results$fullpars[c(9,10,12,13,15,16,18,19)], nrow=4, ncol=2, byrow=T)
  
  # individual parameters
  results$indpar <- data.frame()
  pcount <- 20
  for (i in 1:N){
    results$indpar <- rbind(results$indpar,
                            data.frame(subject=i,
                              kappa=results$fullpars[pcount+1],
                              lambda=results$fullpars[pcount+2],
                              bx1=1,
                              by1=results$fullpars[pcount+3],
                              a1=results$fullpars[pcount+4],
                              bx2=results$fullpars[pcount+5],
                              by2=1,
                              a2=results$fullpars[pcount+6]) )
    pcount <- pcount+6
  }
  
  # make the results object class "grt_wind_fit"
  class(results) <- "grt_wind_fit"
  
  results$N <- N
  results$model <- paste("GRT-wIND_", model, sep="")
  results$predicted <- fitted(results)
  results$observed <- c()
  
  for (i in 1:N){
    results$observed <- c(results$observed, as.vector(pmatrix(cmats[[i]])))    
  }
  
  # get R2 if possible
  tryCatch(
    {
      results$R2 <- cor(results$predicted, results$observed)^2
    },
    error=function(e) {
      results$R2 <- NA
    }
  )
  
  return(results)
  
}




