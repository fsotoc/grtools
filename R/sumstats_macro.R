#' Perform a summary statistics macroanalysis of identification data
#' 
#' Performs a summary statistics (i.e., MSDA) macroanalysis (see Kadlec, 1995;
#' Kadlec & Townsend, 1992) of data from a 2x2 identification experiment. This
#' analysis should be performed together with the microanalysis implemented by
#' the function \code{\link{sumstats_micro}}
#' 
#' @param cmat A 4x4 confusion matrix (see Details).
#' @param use_kadlec If TRUE (default), uses a definition of the decision bound
#'   parameter c from Kadlec (1999). If FALSE, it uses the definition from
#'   MacMillan and Creelman (2005).
#' 
#' @return An object of class "\code{sumstats_macro}"
#' 
#' The function \code{summary} is used to obtain a summary of conclusions from 
#' the analysis about perceptual and decisional separability (see Table 1 in
#' Kadlec, 1995).
#' 
#' @details For an introductory tutorial on the summary statistics macroanalyses, see
#' Ashby & Soto (2005), particularly pages 22-28. 
#' 
#' A 2x2 identification experiment involves two dimensions, A and B,
#' each with two levels, 1 and 2. Stimuli are represented by their level in each
#' dimension (A1B1, A1B2, A2B1, and A2B2) and so are their corresponding correct
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
#' @references
#' Ashby, F. G., & Soto, F. A. (2015). Multidimensional signal detection theory.
#' In J. R. Busemeyer, J. T. Townsend, Z. J. Wang, & A. Eidels (Eds.),
#' \emph{Oxford handbook of computational and mathematical psychology} (pp.
#' 13-34). Oxford University Press: New York, NY.
#' 
#' Kadlec, H. (1995). Multidimensional signal detection analyses (MSDA)
#' for testing separability and independence: A Pascal program. \emph{Behavior
#' Research Methods, Instruments, & Computers, 27}(4), 442-458.
#' 
#' Kadlec, H., & Townsend, J. T. (1992). Signal detection analyses of
#' multidimensional interactions. In F. G. Ashby (Ed.), \emph{Multidimensional
#' models of perception and cognition} (pp. 181–231). Hillsdale, NJ: Erlbaum.
#' 
#' Macmillan, N. A., & Creelman, D. (2005). \emph{Detection theory: A user’s
#' guide (2nd ed.)}. Mahwah, NJ: Erlbaum.
#' 
#' @examples
#' # Create a confusion matrix
#' # Inside the c(...) below, we enter the data from row 1 in the 
#' # matrix, then from row 2, etc.
#' cmat <- matrix(c(140, 36, 34, 40,
#'                  89, 91, 4, 66,
#'                  85, 5, 90, 70,
#'                  20, 59, 8, 163),
#'                  nrow=4, ncol=4, byrow=TRUE)
#'                   
#' # Perform the summary statistics macroanalysis
#' macro_results <- sumstats_macro(cmat)
#' 
#' # See a summary of the results
#' summary(macro_results)
#' 
#' # Print to screen the details of each test
#' macro_results
#' 
#' @seealso \code{\link{sumstats_micro}}
#' 
#' @export
sumstats_macro <- function(cmat, use_kadlec=T) {
  
  ### MACROANALYSES
  macro = list()
  
  # marginal response invariance
  macro$marginal_response_invariance <- mri_test(cmat)
  
  # equal marginal d'
  macro$equal_marginal_d_prime <- emdprime(cmat)
  
  # equal marginal c
  macro$equal_marginal_c <- emc(cmat, use_kadlec=use_kadlec)
 
  # return an object of class "sumstats_macro"
  class(macro) <- "sumstats_macro"
  return(macro)
  
}