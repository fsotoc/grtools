#' Run likelihood ratio tests of separability and independence
#' 
#' Run likelihood ratio tests of perceptual separability, decisional separability,
#' and perceptual independence (see Ashby & Soto, 2015), using the
#' maximum likelihood parameter estimates stored in a \code{grt_wind_fit}
#' object.
#' 
#' @usage lr_test(fitted_model, cmats, n_reps=20, test=c("PS(A)", "PS(B)", "PI", "DS(A)", "DS(B)"))
#' 
#' @param fitted_model An object returned by \code{\link{grt_wind_fit}} or 
#'   \code{\link{grt_wind_fit_parallel}}
#' @param cmats List of confusion matrices. Each entry in the list should contain
#'   the 4x4 confusion matrix from one individual. For a detailed description of
#'   how to create this list, see \code{\link{grt_wind_fit}}
#' @param n_reps  Number of times the optimization algorithm should be run when
#'   a restricted model is fitted to the data (see Details).
#' @param test A string array indicating the test(s) to be performed. "PS(A)" and
#'   "PS(B)" indicate to test perceptual separability of A and B, respectively. "PI"
#'   indicates to test perceptual independence. "DS(A)" and "DS(B)" indicate to test
#'   decisional separability of A and B, respectively. The default is test=c("PS(A)",
#'   "PS(B)", "PI", "DS(A)", "DS(B)"), which includes all tests.
#' 
#' @return An object of class "\code{grt_wind_fit}," including information about
#' likelihood ratio tests.
#' 
#' The function \code{summary} is used to obtain a summary of the model fit to 
#' data and the results of likelihood ratio tests.
#' 
#' @details Each likelihood ratio test involves fitting a restricted model to the
#'   data (e.g., a model in which parameters are fixed to reflect PS) and then 
#'   statistically comparing the fit of the restricted model against that of the
#'   full model (for details, see Ashby & Soto, 2015). The \code{lr_test} function
#'   calls \code{\link{grt_wind_fit_parallel}} to fit a restricted model several
#'   times (determined by parameter \code{n_reps}; set to 20 by default). Each time,
#'   the starting parameter values are the values previously found by fitting the 
#'   full model, with small random values added or subtracted (for details about
#'   the fitting procedure, see \code{\link{grt_wind_fit}}). The best-fitting model
#'   is chosen and used in the likelihood ratio test.
#'   
#'   Because likelihood ratio tests require fitting a GRT-wIND model many times,
#'   you should expect the analysis to take considerable time to finish. We recommend
#'   you to run only the tests that interest you, and not all the tests included by
#'   default.
#' 
#' @references Ashby, F. G., & Soto, F. A. (2015). Multidimensional signal
#'   detection theory. In J. R. Busemeyer, J. T. Townsend, Z. J. Wang, & A.
#'   Eidels (Eds.), \emph{Oxford handbook of computational and mathematical
#'   psychology} (pp. 13-34). Oxford University Press: New York, NY.
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
#' # fit the model to data
#' fitted_model <- grt_wind_fit(cmats)
#' 
#' #' # run the likelihood ratio tests
#' fitted_model <- lr_test(fitted_model, cmats)
#' 
#' # see the results
#' summary(fitted_model)
#' 
#' @export
lr_test <- function(x, ...) UseMethod("lr_test")

#' @export
lr_test.grt_wind_fit <- function(fitted_model, cmats, n_reps=20, test=c("PS(A)", "PS(B)", "PI", "DS(A)", "DS(B)")){
  
  restricted <- list()
  lr_test <- data.frame()
  #---------------------------------------------------------
  # PS(A)
  if (any(test=="PS(A)")) {
    # fit restricted model
    restricted[["PS(A)"]] <- grt_wind_fit_parallel(cmats, start_params = fitted_model$fullpars, model="PS(A)", n_reps=n_reps)
    
    # run likelihood ratio test
    Chi2 <- 2*(-fitted_model$value + restricted[["PS(A)"]]$value)
    dof <- length(fitted_model$par) - length(restricted[["PS(A)"]]$par)
    pval <- pchisq(Chi2, dof, lower.tail = F)
    if (pval<.05){
      viol <- "YES"
    } else {
      viol <- "NO"
    }
    lr_test <- rbind(lr_test, data.frame(Test="Perceptual Separability of dimension A", Chi2=Chi2, DF=dof, 
                                         pval=pval, Violation=viol))
  }
  
  #---------------------------------------------------------
  # PS(B)
  if (any(test=="PS(B)")) {
    # fit restricted model
    restricted[["PS(B)"]] <- grt_wind_fit_parallel(cmats, start_params = fitted_model$fullpars, model="PS(B)", n_reps=n_reps)
    
    # run likelihood ratio test
    Chi2 <- 2*(-fitted_model$value + restricted[["PS(B)"]]$value)
    dof <- length(fitted_model$par) - length(restricted[["PS(B)"]]$par)
    pval <- pchisq(Chi2, dof, lower.tail = F)
    if (pval<.05){
      viol <- "YES"
    } else {
      viol <- "NO"
    }
    lr_test <- rbind(lr_test, data.frame(Test="Perceptual Separability of dimension B", Chi2=Chi2, DF=dof, 
                                         pval=pval, Violation=viol))
  }
  
  
  #---------------------------------------------------------
  # PI
  if (any(test=="PI")) {
    # fit restricted model
    restricted[["PI"]] <- grt_wind_fit_parallel(cmats, start_params = fitted_model$fullpars, model="PI", n_reps=n_reps)
    
    # run likelihood ratio test
    Chi2 <- 2*(-fitted_model$value + restricted[["PI"]]$value)
    dof <- length(fitted_model$par) - length(restricted[["PI"]]$par)
    pval <- pchisq(Chi2, dof, lower.tail = F)
    if (pval<.05){
      viol <- "YES"
    } else {
      viol <- "NO"
    }
    lr_test <- rbind(lr_test, data.frame(Test="Perceptual Independence", Chi2=Chi2, DF=dof, 
                                         pval=pval, Violation=viol))
  }
  
  #---------------------------------------------------------
  # DS(A)
  if (any(test=="DS(A)")) {
    # fit restricted model
    restricted[["DS(A)"]] <- grt_wind_fit_parallel(cmats, start_params = fitted_model$fullpars, model="DS(A)", n_reps=n_reps)
    
    # run likelihood ratio test
    Chi2 <- 2*(-fitted_model$value + restricted[["DS(A)"]]$value)
    dof <- length(fitted_model$par) - length(restricted[["DS(A)"]]$par)
    pval <- pchisq(Chi2, dof, lower.tail = F)
    if (pval<.05){
      viol <- "YES"
    } else {
      viol <- "NO"
    }
    lr_test <- rbind(lr_test, data.frame(Test="Decisional Separability of dimension A", Chi2=Chi2, DF=dof, 
                                         pval=pval, Violation=viol))
  }
  
  #---------------------------------------------------------
  # DS(B)
  if (any(test=="DS(B)")) {
    # fit restricted model
    restricted[["DS(B)"]] <- grt_wind_fit_parallel(cmats, start_params = fitted_model$fullpars, model="DS(B)", n_reps=n_reps)
    
    # run likelihood ratio test
    Chi2 <- 2*(-fitted_model$value + restricted[["DS(B)"]]$value)
    dof <- length(fitted_model$par) - length(restricted[["DS(B)"]]$par)
    pval <- pchisq(Chi2, dof, lower.tail = F)
    if (pval<.05){
      viol <- "YES"
    } else {
      viol <- "NO"
    }
    lr_test <- rbind(lr_test, data.frame(Test="Decisional Separability of dimension B", Chi2=Chi2, DF=dof, 
                                         pval=pval, Violation=viol))
  }
  
  # return restricted models in grt_wind_fit object
  # we keep the results from any previously ran tests
  if (is.null(fitted_model$restricted)){
    fitted_model$restricted <- restricted
    fitted_model$lr_test <- lr_test
  } else {
    fitted_model$restricted <- c(fitted_model$restricted, restricted)
    fitted_model$lr_test <- rbind(fitted_model$lr_test, lr_test)
  }
  
  return(fitted_model)
  
}