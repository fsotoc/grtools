#' Perform a summary statistics microanalysis of identification data
#' 
#' Performs a summary statistics (i.e., MSDA) microanalysis (see Kadlec, 1995;
#' Kadlec & Townsend, 1992) of data from a 2x2 identification experiment. This
#' analysis should be performed together with the macroanalysis implemented by
#' the function \code{\link{sumstats_macro}}
#' 
#' @param cmat A 4x4 confusion matrix (see Details).
#' @param use_kadlec If TRUE (default), uses a definition of the decision bound
#'   parameter c from Kadlec (1999). If FALSE, it uses the definition from
#'   MacMillan and Creelman (2005).
#' 
#' @return An object of class "\code{sumstats_macro}"
#' 
#' The function \code{summary} is used to obtain a summary of conclusions from 
#' the analysis about perceptual independence and decisional separability, (see
#' Table 2 in Kadlec, 1995).
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
#' # Perform the summary statistics microanalysis
#' micro_results <- sumstats_micro(cmat)
#' 
#' # See a summary of the results
#' summary(micro_results)
#' 
#' # Print to screen the details of each test
#' micro_results
#' 
#' @seealso \code{\link{sumstats_macro}}
#' 
#' @export
sumstats_micro <- function(cmat, use_kadlec=T) {
  # sumstats.micro(cmat, use_kadlec=T) produces a GRT summary statistics microanalysis of a 2x2 identification task
  #
  # cmat is a 4x4 confusion matrix with stimuli in rows and responses in columns, each cell has the 
  # frequency of responses for the stimulus/response pair
  #
  # use_kadlec determines whether to use the definition of c proposed by Kadlec (use_kadlec=T, default)
  # or to use the c proposed by MacMillan and Creelman.
  # 
  # Rows and columns should be ordered in the following way: 
  #  	Stimulus/Row 1: A1B1			Response/Column 1: a1b1
  #		Stimulus/Row 2: A2B1			Response/Column 2: a2b1
  #		Stimulus/Row 3: A1B2			Response/Column 3: a1b2
  #		Stimulus/Row 4: A2B2			Respones/Column 4: a2b2
  #
  #
  # The output is a list 
  #
  # Here is an example of the correct way to use this function. At the command window, type the following
  # to create the confusion matrix cmat:
  #
  # > cmat<-matrix(c(140,36,34,40,89,91,4,66,85,5,90,70,20,59,8,163),nrow=4,ncol=4,byrow=TRUE)
  #
  # Then you can call the function by typing:
  #
  # > sumstats.micro(cmat)
  #
  #
  # Fabian A. Soto
  # Department of Psychological and Brain Sciences
  # University of California, Santa Barbara
  # February 2, 2015
  

  
  micro <- list()
  
  # sampling independence
  micro$sampling_independence <- sind(cmat)
  
  # equal conditional d's
  micro$equal_conditional_d_prime <- econdprime(cmat)
  
  # equal conditional c
  micro$equal_conditional_c <- econdc(cmat)
  

  # return object of class "sumstats_micro"
  class(micro) <- "sumstats_micro"
  return(micro)
}