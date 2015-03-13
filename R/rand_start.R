# adds a random value to starting parameters and checks that it is within
# acceptable limits
rand_start<-function(start_params, low_params, high_params,rand_pert){
  cand_params=runif(length(start_params),-1,1)*rand_pert+start_params
  while(any(cand_params>high_params) || any(cand_params<low_params)){
    cand_params=runif(length(start_params),-1,1)*rand_pert+start_params
  }
  return (cand_params)
}