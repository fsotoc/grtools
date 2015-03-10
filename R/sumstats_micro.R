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