matrix_sample<-function(p_matrix, N_row){
  # samples a confusion matrix (frequency of responses) from a matrix of confusion 
  # probabilities using the multinomial distribution
  # prob_matrix is a matrix of confusion probabilities
  # N_row is the number of trials per row; it can be a single number (equal number of
  # trials per row) or a vector with four numbers (unequal number of trials per row)
  
  if ( length(N_row)==1 ){
    N_row = rep(N_row, 4)
  } else if ( length(N_row)!=4 ){
    stop("N_row should be size 1 or 4")
  }
  
  f_matrix <- matrix(0, nrow=4, ncol=4, byrow=T)
  for (i in 1:4){
    f_matrix[i,] <- rmultinom(1, N_row[i], p_matrix[i,])
  }
  
  return(f_matrix)
  
}