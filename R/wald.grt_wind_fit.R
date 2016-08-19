#' Run Wald tests of separability and independence
#' 
#' Run Wald tests of perceptual separability, decisional separability, and
#' perceptual independence (see appendix of Soto et al., 2015), using the
#' maximum likelihood parameter estimates stored in a \code{grt_wind_fit}
#' object.
#' 
#' @usage wald(fitted_model, data, estimate_hess=T)
#' 
#' @param fitted_model An object returned by \code{\link{grt_wind_fit}}
#' @param data List of confusion matrices. Each entry in the list should contain
#'   the 4x4 confusion matrix from one individual. For a detailed description of
#'   how to create this list, see \code{\link{grt_wind_fit}}
#' @param estimate_hess Indicates whether a numerical estimate of the Hessian
#'   should be obtained. If an estimate has been previously obtained and is
#'   stored in \code{fitted_model}, setting \code{estimate_hess=F} will save
#'   computing time.
#' 
#' @return An object of class "\code{grt_wind_fit}," including information about
#' Wald tests.
#' 
#' The function \code{summary} is used to obtain a summary of the model fit to 
#' data and the results of Wald tests.
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
#' data <- list(matrix(c(100,1,9,8,10,110,7,4,31,3,80,10,54,4,52,19),nrow=4,ncol=4,byrow=TRUE))
#' data[[2]] <- matrix(c(122,7,0,1,1,102,1,5,3,2,111,9),nrow=4,ncol=4,byrow=TRUE)
#' data[[3]] <- matrix(c(107,0,5,0,0,101,1,4,3,1,113,2,0,3,1,108),nrow=4,ncol=4,byrow=TRUE)
#' data[[4]] <- matrix(c(122,1,0,0,1,120,0,0,0,0,118,6,0,1,6,118),nrow=4,ncol=4,byrow=TRUE)
#' data[[5]] <- matrix(c(89,17,6,4,4,81,8,6,14,7,86,1,11,25,17,26),nrow=4,ncol=4,byrow=TRUE)
#' 
#' # fit the model to data
#' fitted_model <- grt_wind_fit(data)
#' 
#' #' # run the wald tests
#' fitted_model <- wald(fitted_model, data)
#' 
#' # see the results
#' summary(fitted_model)
#' 
#' @export
wald <- function(x, ...) UseMethod("wald")

#' @export
wald.grt_wind_fit <- function(fitted_model, data, estimate_hess=T){
  
  if (fitted_model$model=="GRT-wIND"){
    # get number of subjects N and parameters K
    N <- length(data)
    K <- length(fitted_model$par)
    
    # get hessian
    if (estimate_hess==T){
      H <- numDeriv::hessian(func=grt_wind_nll, x=fitted_model$par, data=data, method.args=list(d=0.01, r=6))
    } else {
      H <- fitted_model$hessian
    }
    
    # estimate covariance matrix by inverting the Hessian
    S <- solve(H)
    if (any(is.infinite(S))){
      S <- MASS::ginv(H)          # if inversion failed, get P-M pseudoinverse (requires package MASS)
    }
    
    #----------------------------------------------
    # global test for PS of dimension A from dimension B
    R <- matrix(0, nrow=4, ncol=K)
    R[1,3] <- 1
    R[2,11] <- 1
    R[3,1] <- 1
    R[3,5] <- -1
    R[4,8] <- 1
    R[4,14] <- -1
    
    q <- matrix(c(0, 1, 0, 0), ncol=1, nrow=4)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- data.frame(Test="Perceptual Separability of dimension A", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol)
    
    
    #---------------------------------------------
    # global test for PS of dimension B from dimension A
    R <- matrix(0,nrow=4,ncol=K)
    R[1,2] <- 1
    R[2,9] <- 1
    R[3,4] <- 1
    R[3,6] <- -1
    R[4,12] <- 1
    R[4,15] <- -1
    
    q <- matrix(c(0, 1, 0, 0), ncol=1, nrow=4)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Perceptual Separability of dimension B", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    #------------------------------------------------
    # global test for PI
    R <- matrix(0, nrow=4, ncol=K)
    R[1,7] <- 1
    R[2,10] <- 1
    R[3,13] <- 1
    R[4,16] <- 1
    
    q <- matrix(c(0, 0, 0, 0), ncol=1, nrow=4)
    wtr <- wald_test(R, q, fitted_model$par, S)  
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Perceptual Independence", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    #----------------------------------------------
    # Test of average DS of dimension A from dimension B
    R <- matrix(0, nrow=1, ncol=K)
    R[16+((1:N)-1)*6+3] <- 1
    
    q <- matrix(0, ncol=1, nrow=1)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Average Decisional Separability of dimension A", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    #------------------------------------------------
    # Test of individual DS of dimension A from dimension B
    R <- matrix(0, nrow=N, ncol=K)
    for (i in 1:N){
      R[i,16+(i-1)*6+3] <- 1
    }
    
    q <- matrix(0, nrow=N, ncol=1)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Individual Decisional Separability of dimension A", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    #----------------------------------------------
    # Test of average DS of dimension B from dimension A
    R <- matrix(0, nrow=1, ncol=K)
    R[16+((1:N)-1)*6+5] <- 1
    
    q <- matrix(0, ncol=1, nrow=1)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Average Decisional Separability of dimension B", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    #------------------------------------------------
    # Test of individual DS of dimension B from dimension A
    R <- matrix(0, nrow=N, ncol=K)
    for (i in 1:N){
      R[i,16+(i-1)*6+5] <- 1
    }
    
    q <- matrix(0, nrow=N, ncol=1)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Individual Decisional Separability of dimension B", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    
    #--------------------------------------------------
    # Test of equal variances
    R <- matrix(0, nrow=6, ncol=K)
    R[1,8] <- 1
    R[2,9] <- 1
    R[3,11] <- 1
    R[4,12] <- 1
    R[5,14] <- 1
    R[6,15] <- 1
    
    q <- matrix(1, nrow=6, ncol=1)
    wtr <- wald_test(R, q, fitted_model$par, S)
    # put results in table
    wt_table <- rbind(wt_table, data.frame(Test="Equal Variances in All Distributions", Wald_Statistic=wtr$stat, df=wtr$df, pvalue=wtr$pval, Violation=wtr$viol))
    
    #---------------------------------------------------
    # change column names
    names(wt_table) <- c("Test", "Wald stat", "DF", "P-value", "Violation?")
    
    # return a new grt_wind_fit object, but with hessian updated
    # and a new results table
    fitted_model$wald_test <- wt_table
    fitted_model$hessian <- H
    
  } else if (fitted_model$model=="GRT-wIND: VAR1"){
    fitted_model <- wald_grt_wind_VAR1(fitted_model, data, estimate_hess)
  }
  
  return(fitted_model)
  
}


