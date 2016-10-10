##' 
##' hessfix. Find hessian with fixed parameters
##'
##' specify an argument 'fixed', a vector of numerical values.
##' created just as in function optifix
##' 
##' @export
##' 

hessfix <- function(par, fixed, fn, ..., method.args=list()) {
  
  .npar=length(par)
  indx = 1:.npar
  .fixValues = par[indx==fixed]
  .fixIndex = fixed&(indx!=fixed)
  .par = par[!fixed]
  
  
  .fn <- function(par,...){
    .par = rep(NA,length(fixed))
    .par[!fixed] = par
    .par[indx==fixed] = .fixValues
    .par[.fixIndex] = .par[fixed[.fixIndex]]
    fn(.par,...)
  }
  
  H = numDeriv::hessian(func=.fn, x=.par, ..., method.args=method.args)
  
  return(H)
}
