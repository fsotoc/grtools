#' Plot a \code{grt_hm_fit} object
#' 
#' Plot the object returned by \code{\link{grt_hm_fit}}
#' 
#' @param model A \code{grt_hm_fit} object
#' @param labels Optional names for the labels of dimensions A and B
#' @param ellipse_width Parameter controlling the width of the drawn ellipses
#' @export
plot.grt_hm_fit <- function(x, labels=c("dim A", "dim B"), ellipse_width=0.8, ...){
  # ellipse_width determines the width of the ellipses
  # labels determines the labels for each axis
  
  model <- x$best_model
  
  
  plot.new()
  
  # first plot the main panel
  par(mar=c(4,4,2,2)+0.1,fig=c(0.2,1,0.2,1), lty=1, new=TRUE)
  
  # get range of values
  ranx <- c(min(model$means[,1]-2), max(model$means,1)+2)
  rany <- c(min(model$means[,2]-2), max(model$means,2)+2)
  
  # draw model$means of distributions
  plot(model$means[,1], model$means[,2], pch=3, 
       xlim=ranx,
       ylim=rany,
       xlab=labels[1],
       ylab=labels[2])
  
  
  # draw contours of distributions
  ellipse <- function(s,t) {u<-c(s,t)-center; u %*% sigma.inv %*% u / 2}
  n <- 200
  x<-1:200/10-10
  y<-1:200/10-10
  for (i in 1:4){
    center <- model$means[i,]
    sigma.inv <- solve(model$covmat[[i]])
    z <- mapply(ellipse, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, `+`)))
    contour(x,y,matrix(z,n,n), levels=ellipse_width,drawlabels=F, add=T)
  }
  
  # add decision bounds
#   } else { 
#     abline(a=model$a1, b=model$by1)
#     abline(a=model$a1, b=model$by1)
#   }
  
    
  # add marginal distributions at the bottom
  par(mar=c(1,3.7,1,1.7)+0.1,fig=c(0.2,1,0,0.2), new=T)
  for (i in 1:4){
    x <- 1:100*(ranx[2]-ranx[1])/100+ranx[1]
    y <- dnorm(x,mean=model$means[i,1],sd=sqrt(model$covmat[[i]][1,1]))
    if (i>2){par(new=T,lty=2)} else {par(new=T,lty=1)}
    plot(x,y,type="l", axes=F, ylab="", xlab="", xlim=ranx)
  }
  par(new=T)
  Axis(side=1)
  
  # add marginal distributions to the left
  par(mar=c(3.7,1,1.7,1)+0.1, fig=c(0,0.2,0.2,1), new=T)
  for (i in 1:4){
    x <- 1:100*(rany[2]-rany[1])/100+rany[1]
    y <- dnorm(x,mean=model$means[i,2],sd=sqrt(model$covmat[[i]][2,2]))
    if (i==2 | i==4){par(new=T,lty=2)} else {par(new=T,lty=1)}
    plot(y,x,type="l", axes=F, ylab="", xlab="", ylim=rany)
  }
  par(new=T)
  Axis(side=2)
  
  # add scatterplot if there are predicted and observed values
  if (any(names(model)=="predicted") & any(names(model)=="observed")) {
    par(mar=c(1.5,1.5,1,1), fig=c(0,0.33,0,0.33), new=T)
    plot(model$predicted,model$observed,pch=21,cex=.3,col='gray40',bg='gray40',bty='n', axes=F)
    abline(a=0, b=1, lty=1)
    axis(side=1, at=c(0,1), mgp=c(3,0.5,0))
    axis(side=2, at=c(0,1), mgp=c(3,0.5,0))
  }
  
}

  
  
