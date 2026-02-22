#' Fit a GRT-wIND model to identification data with 3 "A" responses and 2 "B" responses
#' 
#' See help file for \code{\link{grt_wind_fit}}
#' 
#' @export
grt_wind_fit_r3x2 <- function(cmats, start_params=c(), model="1-VAR", rand_pert=0.3, control=list()) {  
  
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
    bx1_start <- rep(0,N)                              # start with bounds orthogonal to dimension (DS)
    bx2_start <- rep(0,N)                              # start with bounds orthogonal to dimension (DS)
    by_start <- rep(0,N)
    ax1_start <- rep(-0.25,N)                           # start with bound right between the two means
    ax2_start <- rep(-0.75,N)                           # start with bound right between the two means
    ay_start <- rep(-0.5,N)
    
    for (i in 1:N) {
      start_params <- c(start_params, kappa_start[i], lambda_start[i], bx1_start[i], ax1_start[i], bx2_start[i], ax2_start[i], by_start[i], ay_start[i])
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
    low_params <- c(low_params, kappa_low, lambda_low, b_low, a_low, b_low, a_low, b_low, a_low)
    up_params <- c(up_params, kappa_up, lambda_up, b_up, a_up, b_up, a_up, b_up, a_up)
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
  
  # GRTapas
  if (model=="GRTapas"){
    # fix distribution of A2B1 to be same as A1B1
    fixed[3] <- 1
    fixed[4] <- 2
    fixed[12] <- 9
    fixed[13] <- 10
    fixed[14] <- 11
    
    # make sure to set starting parameters to values that make sense
    # means with mickey mouse pattern
    start_params[1:8] <- c(0,0,0,0,-1,1,1,1)
    
    # means cannot be negative in this model along the y-axis
    low_params[c(2,4,6,8)] <- 0
    # deal with individual parameters
    for (i in 1:N) {

      # we start bounds in dimension A at -2 and 2
      start_params[20+((i-1)*8)+4] <- 2
      start_params[20+((i-1)*8)+6] <- -2
      
    }
    
  }
  
  # PS(A)
  if (any(model=="PS(A)")){
    # fix means along x axis
    fixed[5] <- 1
    fixed[7] <- 3
    # fix variances along x axis
    fixed[15] <- 9
    fixed[18] <- 12
  }
  # PS(B)
  if (any(model=="PS(B)")){
    # fix means along y axis
    fixed[4] <- 2
    fixed[8] <- 6
    # fix variances along y axis
    fixed[13] <- 10
    fixed[19] <- 16
  } 
  # PI
  if (any(model=="PI")) {
    # make sure that the starting value of correlations is zero
    start_params[c(11, 14, 17, 20)] <- 0
    # set value of correlations to starting value
    fixed[c(11, 14, 17, 20)] <- c(11, 14, 17, 20)
  } 
  # DS(A)
  if (any(model=="DS(A)")){
    for (i in 1:N){
      # make sure that the starting bounds are orthogonal
      start_params[20+(i-1)*8+3] <- 0
      start_params[20+(i-1)*8+5] <- 0
      # set value of the bound slope to starting value
      fixed[20+(i-1)*8+3] <- 20+(i-1)*8+3
      fixed[20+(i-1)*8+5] <- 20+(i-1)*8+5
    }
  } 
  # DS(B)
  if (any(model=="DS(B)")){
    for (i in 1:N){
      # make sure that the starting bound is orthogonal
      start_params[20+(i-1)*8+7] <- 0
      # set value of the bound slope to starting value
      fixed[20+(i-1)*8+7] <- 20+(i-1)*8+7
    }
  } 
  # 1-VAR
  if (any(model=="1-VAR")){
    # make sure that all starting variances are equal to one
    start_params[c(9,10,12,13,15,16,18,19)] <- 1
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
  results <- optifix(start_params, fixed, grt_wind_nll_r3x2, data=cmats, method='L-BFGS-B', lower=low_params, upper=up_params, control=control, hessian = F)
  
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
  # parameters for bounds of x-axis are indexed 1 and 2 in names
  # parameters for single bound of y-axis are indexed 3 in names
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
                              bx2=1,
                              by2=results$fullpars[pcount+5],
                              a2=results$fullpars[pcount+6],
                              bx3=results$fullpars[pcount+7],
                              by3=1,
                              a3=results$fullpars[pcount+8]) )
    pcount <- pcount+8
  }
  
  # make the results object class "grt_wind_fit"
  class(results) <- "grt_wind_fit"
  
  results$N <- N
  results$model <- "GRT-wIND"
  results$restrictions <- c(model, "r3x2")
  results$predicted <- fitted_grt_wind_fit_r3x2(results)
  results$observed <- c()
  
  for (i in 1:N){
    results$observed <- c(results$observed, as.vector(pmatrix(cmats[[i]])))    
  }
  
  # get R2 if possible
  tryCatch(
    {
      results$R2 <- cor(results$predicted[which(!is.na(results$observed))], 
                        results$observed[which(!is.na(results$observed))])^2
    },
    error=function(e) {
      results$R2 <- NA
    }
  )
  
  return(results)
  
}




