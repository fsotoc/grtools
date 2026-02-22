#' Calculate confidence intervals from a matrix of parameter samples
#'
#' Takes as input a NxK matrix of parameter samples, with K parameters in columns and N 
#' samples in rows and, for each column, it computes a (1-alpha)\% two-sided confidence interval 
#' for that parameter using a percentile method.
#' 
#' @param input_matrix The NxK input matrix, composed of N rows of K parameter vectors.
#' @param mle_pars Optional vector of parameter estimates associated with each
#' selected column. If omitted, estimates are computed from \code{input_matrix}
#' as column means (or the scalar mean for a vector input).
#' @param col_selection An optional vector of column indexes (e.g., \code{c(1,2,3,7:20)})
#' or column names (e.g., \code{c("A1B1_mean_A","A2B1_mean_B")}) for which the confidence interval
#' will be computed. If not provided, confidence intervals for all parameters are reported.
#' @param alpha Error probability parameter used to create (1-alpha)*100\% two-sided confidence
#' interval. Defaults to 0.05 (95\% confidence interval).
#' @return A matrix including maximum likelihood estimates of each selected parameter, plus
#' lower and upper limits of the confidence interval.
#' 
#' @details The percentile method is used to calculate confidence intervals from a set of samples.
#' For more details, see Carpenter and Bithell (2000).
#' 
#' @references 
#' Carpenter, J., & Bithell, J. (2000). Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians. Statistics in medicine, 19(9), 1141-1164.
#' 
#' @examples
#' # create a matrix of data that is ordered columnwise:
#' data <- matrix(c(161,41,66,32,24,147,64,65,37,41,179,43,14,62,22,202,
#' 26,82,67,25,8,188,54,50,34,75,172,19,7,103,14,176,
#' 117,64,89,30,11,186,69,34,21,81,176,22,10,98,30,162,
#' 168,57,47,28,15,203,33,49,58,54,156,32,9,96,9,186,
#' 169,53,53,25,34,168,69,29,38,48,180,34,19,44,60,177),nrow=5, ncol=16, byrow=TRUE)
#' 
#' #calculating 95\% confidence interval for all data
#' ci_data <- confidence_interval(data, alpha=0.05)
#' 
#' 
#' @export
confidence_interval <- function(input_matrix, mle_pars = NULL, col_selection = NULL, alpha = 0.05){
  cn <- NULL
  
  if (is.null(dim(input_matrix))) {
    input_matrix <- matrix(input_matrix, ncol = 1)
  }
  
  # if no columns were selected, do all of them
  if (is.null(col_selection)){
    col_selection <- 1:ncol(input_matrix)
  }
  
  # keep column names if given
  if (is.character(colnames(input_matrix))){
    if (is.character(col_selection)){
      cn <- col_selection
    } else {
      cn <- colnames(input_matrix)[col_selection]
    }
  }
  
  input_matrix <- input_matrix[, col_selection, drop = FALSE]
  
  # deal with case of a single vector
  if (is.vector(input_matrix) || (is.matrix(input_matrix) && ncol(input_matrix) == 1)) {
    vec <- as.vector(input_matrix)
    # get parameter estimates and bounds of the confidence interval
    mle <- if (is.null(mle_pars)) mean(vec) else mle_pars
    ci_low <- quantile(vec, probs=alpha/2)
    ci_high <- quantile(vec, probs=1-alpha/2)
    if (is.null(cn)) {
      cn <- if (!is.null(colnames(input_matrix))) colnames(input_matrix) else "parameter"
    }
  } else {
    # get parameter estimates and bounds of the confidence interval
    mle <- if (is.null(mle_pars)) colMeans(input_matrix) else mle_pars
    ci_low <- apply(input_matrix, 2, quantile, probs=alpha/2)
    ci_high <- apply(input_matrix, 2, quantile, probs=1-alpha/2)
    if (is.null(cn)) {
      cn <- colnames(input_matrix)
    }
  }
  
  final_matrix <- matrix(c(mle, ci_low, ci_high), nrow = 3, byrow=T)
  colnames(final_matrix) <- cn
  rownames(final_matrix) <- c("Estimate:", 
                              paste((1-alpha)*100,"% confidence interval, lower bound:", sep=""), 
                              paste((1-alpha)*100,"% confidence interval, upper bound:", sep=""))
  
  return(final_matrix)
}
