#' @export
data_sample <- function(x, ...) UseMethod("data_sample")

#' @export
data_sample.grt_wind_fit <- function(model, N_row){
  # samples a list of confusion matrices from a grt-wind model
  # model is a list either obtained directly from grt_wind_fit, or created manually
  # and declared of class "grt_wind_fit". See the function grt_wind_fit to look
  # for the required list components
  # N_row is the number of trials per row; it can be a single number (equal number of
  # trials per row; same for all participants), a vector with four numbers (unequal number 
  # of trials per row, but same across participants), or a list (one element per 
  # participant) of such single values or four-element vectors.
  
  
  # get a list of confusion probability matrices from the model
  pmats <- fitted(model, out_format="list")
  
  # get a list of confusion matrices (response frequencies)
  cmats <- list()
  for (sub in 1:model$N){
    # if a list of N_rows was provided, get value for this participant
    if (class(N_row)=="list"){
      nr <- N_row[[sub]]
    } else {
      nr <- N_row
    }
    
    cmats <- c(cmats, list(matrix_sample(pmats[[sub]], nr)))
  }
  
  return(cmats)
  
}