##' 
##' optifix. Optimise with fixed parameters
##'
##' its like optim, but with fixed parameters.
##' 
##' specify a second argument 'fixed', a vector of numerical values.
##' There are three options
##' (a) To indicate that the parameter in fn() is optimised over, use 0
##' (b) To indicate that the parameter in fn() is fixed to its starting value, 
##' use its index. Thus, if the Nth parameter in fn() is fixed to the starting
##' value, then fixed[N] <- N.
##' (c) To indicate that the parameter in fn() is fixed to the value of another
##' parameter, use the index of that parameter. Thus, if the Nth parameter in fn()
##' should have the same value as the Mth parameter, then fixed[N] <- M.
##'
##' The return thing is the return thing from optim() but with a couple of extra
##' bits - a vector of all the parameters and a vector copy of the 'fixed' argument.
##'
##' This function was adapted by Fabian Soto from code originally written by 
##' Barry Rowlingson <b.rowlingson@lancaster.ac.uk> in October 2011 and released 
##' under a CC By-SA license: http://creativecommons.org/licenses/by-sa/3.0/
##'
##' Any use of this script must retain the text: 
##' "Originally written by Barry Rowlingson" in comments.
##' 
##' @export
optifix <- function(par, fixed, fn, gr = NULL, ...,
           method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
           lower = -Inf, upper = Inf,
           control = list(), hessian = FALSE){
   force(fn)
   force(fixed)
   .npar=length(par)
   indx = 1:.npar
   .fixValues = par[indx==fixed]
   .fixIndex = fixed&(indx!=fixed)

   .parStart = par[!fixed]
   .lower = lower[!fixed]
   .upper = upper[!fixed]
   
   .fn <- function(par,...){
     .par = rep(NA,length(fixed))
     .par[!fixed] = par
     .par[indx==fixed] = .fixValues
     .par[.fixIndex] = .par[fixed[.fixIndex]]
     fn(.par,...)
   }


# usually we don't provide a gradient, so this is commented out

#    if(!is.null(gr)){
#      .gr <- function(par,...){
#        .gpar = rep(NA,sum(!fixed))
#        .gpar[!fixed] = par
#        .gpar[fixed] = .fixValues
#        gr(.gpar,...)[!fixed]
#      }
#    }else{
     .gr <- NULL
   # }

   .opt = optim(.parStart,.fn, .gr,..., method=method, lower=.lower, upper=.upper,
                control=control, hessian=hessian)

   .opt$fullpars = rep(NA,.npar)
   .opt$fullpars[indx==fixed]=.fixValues
   .opt$fullpars[!fixed]=.opt$par
   .opt$fullpars[.fixIndex]=.opt$fullpars[fixed[.fixIndex]]
   .opt$fixed = fixed
   return(.opt)
   
 }
