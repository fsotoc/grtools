pmatrix <- function(m) {
  
  nrows <- dim(m)[1]
  ncols <- dim(m)[2]
  
  #define a matrix to be filled with proportions
  y <- matrix(data=NA, ncol=ncols, nrow=nrows, byrow=T)
  
  # computing proportions, by row, assign proportions to matrix
  for (i in 1:nrows){
    y[i,] <- m[i,] / sum(m[i,])  
  }

  #output
  return(y)
}