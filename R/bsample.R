# function that performs parametric bootstrap sampling from a fitted grt_wind model
# see bootstrap_sample for help about inputs and outputs

bsample<-function(fitted_model, N_rows, N){  
  
  # preallocating a matrix with proper dimensions to store samples
  # rows are samples and columns are model parameters
  BS <- matrix(nrow = N, ncol=length(fitted_model$fullpars), byrow = TRUE)
  
  # for each sample
  for (i in 1:N){
    
    # sample data from the original fitted GRT-wIND model
    sampled_cmats <- data_sample(fitted_model, N_rows)
    
    # fit model to sampled data
    # we start the optimization procedure at the maximum likelihood estimates without any
    # random perturbation
    bootstrap_model <- grt_wind_fit(sampled_cmats, start_params = fitted_model$fullpars, 
                                    model=fitted_model$restrictions, rand_pert = 0)  
    
    # place estimated parameters in matrix with bootstrap samples
    BS[i,] <- bootstrap_model$fullpars
  }
  
  # names of the parameters for the group model
  paramether1 <- c("A1B1_mean_A", "A1B1_mean_B", "A2B1_mean_A", "A2B1_mean_B", 
                   "A1B2_mean_A", "A1B2_mean_B", "A2B2_mean_A", "A2B2_mean_B",
                   "A1B1_var_A", "A1B1_var_B", "A1B1_corr","A2B1_var_A", "A2B1_var_B", "A2B1_corr",
                   "A1B2_var_A", "A1B2_var_B", "A1B2_corr","A2B2_var_A", "A2B2_var_B", "A2B2_corr")
  
  # names of the individual model parameters
  paramether2 <- matrix(vector (mode="character",length=NCOL(BS)-20), ncol = 1, byrow = TRUE)
  J<-c()
  # parameters defined for each subject
  J[1]="_kappa"
  J[2]="_lambda"
  J[3]="_by1"
  J[4]="_a1"
  J[5]="_bx2"
  J[6]="_a2"
  J[7]="_by2"
  J[8]="_bx3"
  J[9]="_a3"
  # number of individual parameters is different for r3x2 model (type-2 GRTapas)
  if (any(fitted_model$restrictions=="r3x2")){
    paramether2 <- paste0("subject",floor((row(paramether2)-1)/8)+1, J[rep(c(1:4,7,6,8,9), fitted_model$N)])
  } else {
    paramether2 <- paste0("subject",floor((row(paramether2)-1)/6)+1, J[(row(paramether2)-1)%%6+1])
  }
  paramether <- append(paramether1, paramether2)
  colnames(BS) <- paramether
  
  return(BS) #returning the final matrix as the output
}