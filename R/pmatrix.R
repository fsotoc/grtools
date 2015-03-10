pmatrix <- function(m) {
  
  #define a matrix to be filled with proportions
  y <- matrix(data=NA, ncol=4, nrow=4, byrow=T)
  
  # computing proportions, by row, assign proportions to matrix
  for (i in 1:4){
    y[i,] <- m[i,] / sum(m[i,])  
  }

  #output
  return(y)
}