#' @export
sumstats_macro <- function(cmat, use_kadlec=T) {
  # sumstats.macro(cmat, use_kadlec=T) produces a GRT summary statistics macroanalysis of a 2x2 identification task
  #
  # cmat is a 4x4 confusion matrix with stimuli in rows and responses in columns, each cell has the 
  # frequency of responses for the stimulus/response pair
  #
  # use_kadlec determines whether to use the definition of c proposed by Kadlec (use_kadlec=T, default)
  # or to use the c proposed by MacMillan and Creelman.
  # 
  # Rows and columns should be ordered in the following way: 
  #		Stimulus/Row 1: A1B1			Response/Column 1: a1b1
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
  # > sumstats.macro(cmat)
  #
  #
  # Fabian A. Soto
  # Department of Psychological and Brain Sciences
  # University of California, Santa Barbara
  # January 29, 2015

  
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