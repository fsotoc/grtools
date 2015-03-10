wald_test <- function(R, q, w, covmat){
  # perform a Wald test defined by the restrictions matrix R,
  # the vector of values for restrictions q, the vector of
  # parameters w and their covariance matrix covmat
  
  dv_covmat <- R%*%covmat%*%t(R)
  dv_infmat <- solve(dv_covmat)
  if (any(is.infinite(dv_infmat))){
    dv_infmat <- MASS:::ginv(dv_covmat)  # if inversion failed, use P-M pseudoinverse (requires package MASS)
  }
  
  
  stat <- t(R%*%w-q)%*%dv_infmat%*%(R%*%w-q)
  df <- dim(R)[1]
  pval <- pchisq(stat, df, lower.tail=F)
  
  if (pval<.05){
    viol <- "YES"
  } else {
    viol <- "NO"
  }
  
  return(data.frame(stat, df, pval, viol))

  
}