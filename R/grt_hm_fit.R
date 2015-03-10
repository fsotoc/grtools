
rand_start<-function(start_params, low_params, high_params,rand_pert){
    cand_params=runif(length(start_params),-1,1)*rand_pert+start_params
    while(any(cand_params>high_params) || any(cand_params<low_params)){
        cand_params=runif(length(start_params),-1,1)*rand_pert+start_params
    }
    return (cand_params)
}

optims<-function(data, rand_pert=0, reps=1){     

    mle_model1=list(value=Inf)
    mle_model2=list(value=Inf)
    mle_model3=list(value=Inf)
    mle_model4=list(value=Inf)
    mle_model5=list(value=Inf)
    mle_model6=list(value=Inf)
    mle_model7=list(value=Inf)
    mle_model8=list(value=Inf)
    mle_model9=list(value=Inf)
    mle_model10=list(value=Inf)
    mle_model11=list(value=Inf)
    mle_model12=list(value=Inf)
    
    start_params=c(1,1,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf)
    for(i in 1:reps){
       maybe=rand_start(start_params,low_params,high_params,rand_pert)
       maybe_mod=optim(par=maybe, fn=negloglik_mod1, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
       if(maybe_mod$value<mle_model1$value){
           start_params<-maybe
       }
       mle_model1<-optim(par=start_params, fn=negloglik_mod1, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,0,1,1,-0.5,-0.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod2, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model2$value){
            start_params<-maybe
        }
        mle_model2<-optim(par=start_params, fn=negloglik_mod2, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 

    
    start_params=c(1,0,1,1,-0.5,-0.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod3, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model3$value){
            start_params<-maybe
        }
        mle_model3<-optim(par=start_params, fn=negloglik_mod3, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,1,0,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod4, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model4$value){
            start_params<-maybe
        }
        mle_model4<-optim(par=start_params, fn=negloglik_mod4, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 

    start_params=c(1,0,1,1,0,-0.5,-0.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod5, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model5$value){
            start_params<-maybe
        }
        mle_model5<-optim(par=start_params, fn=negloglik_mod5, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
        
    start_params=c(1,0,0,1,1,1,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod6, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model6$value){
            start_params<-maybe
        }
        mle_model6<-optim(par=start_params, fn=negloglik_mod6, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,0,1,1,0,-0.5,-0.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod7, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model7$value){
            start_params<-maybe
        }
        mle_model7<-optim(par=start_params, fn=negloglik_mod7, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,1,0,0,0,0,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,1,1,1,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod8, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model8$value){
            start_params<-maybe
        }
        mle_model8<-optim(par=start_params, fn=negloglik_mod8, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,0,1,1,0,0,0,0,-0.5,-0.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod9, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model9$value){
            start_params<-maybe
        }
        mle_model9<-optim(par=start_params, fn=negloglik_mod9, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,0,0,1,1,1,0,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod10, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model10$value){
            start_params<-maybe
        }
        mle_model10<-optim(par=start_params, fn=negloglik_mod10, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,0,1,1,0,0,0,0,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf)
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod11, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model11$value){
            start_params<-maybe
        }
        mle_model11<-optim(par=start_params, fn=negloglik_mod11, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
    
    start_params=c(1,0,0,1,1,1,0,0,0,0,-.5,-.5)
    grad_step <- rep(1e-1, times=length(start_params))
    low_params=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
    high_params=c(Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf) 
    for(i in 1:reps){
        maybe=rand_start(start_params,low_params,high_params,rand_pert)
        maybe_mod=optim(par=maybe, fn=negloglik_mod12, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
        if(maybe_mod$value<mle_model12$value){
            start_params<-maybe
        }
        mle_model12<-optim(par=start_params, fn=negloglik_mod12, data=data, method="L-BFGS-B", lower=low_params, upper=high_params, control=list(ndeps=grad_step))
    } 
        
    optim_list<-list(mle_model1,mle_model2,mle_model3,mle_model4,mle_model5,mle_model6,mle_model7,mle_model8,mle_model9,mle_model10,mle_model11,mle_model12)
    return(optim_list)
}

order_aic <-function(optim_list) {
    aic_list=rep(0,12)
    L=rep(0,12)
    model=c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", "{PI, DS}", "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", "{1_RHO, DS}", "{PS(B), DS}", "{DS}")
    for(i in 1:12){
        L[i]=-optim_list[[i]]$value
        m=length(optim_list[[i]]$par)
       aic_list[i] = -2*L[i]+2*m+(2*m^2+2*m)/(16-m-1)
    }

    aic_exp=rep(0,12)    
    aic_weight=rep(0,12)
    aic_exp <- exp(-(aic_list-min(aic_list))/2)
    aic_weight <- aic_exp/sum(aic_exp)
    
ordered_aic=data.frame(model,L,aic_list,aic_weight)
colnames(ordered_aic)=c("model","log-likelihood", "AIC", "AIC weight")
ordered_aic=ordered_aic[order(-aic_weight),]
ordered_aic[4]=prettyNum(round(ordered_aic[4],digits=3),nsmall=2)
return(ordered_aic)
}

#' @export
grt_hm_fit<-function(data, rand_pert=0, reps=1){
    optim_list=optims(data,rand_pert,reps)
    o_aic=order_aic(optim_list)
    
    best_model=list()
    best_model$means<-matrix(0,4,2,byrow=TRUE)
    best_model$covmat=list()
    best_model$a1=0
    best_model$a2=0
    best_model$model=""
    
    o_match=o_aic[1,]$model
    best_model$model=paste("GRT-", o_match,sep="")
    model_list=c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", "{PI, DS}", "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", "{1_RHO, DS}", "{PS(B), DS}", "{DS}")
    model_num=pmatch(o_match,model_list)
    w=optim_list[[model_num]]$par
    best_model$model=paste("GRT-", o_match,sep="")
    switch(model_num,
           mod1={best_model$means[2,1]<-w[1]
                 best_model$means[3,2]<-w[2]
                 best_model$means[4,1]<-w[1]
                 best_model$means[4,2]<-w[2]
                 best_model$covmat[[1]]<-diag(2)
                 best_model$covmat[[2]]<-diag(2)
                 best_model$covmat[[3]]<-diag(2)
                 best_model$covmat[[4]]<-diag(2)
                 best_model$a1=w[3]
                 best_model$a2=w[4]},
           
           mod2={best_model$means[2,1]<-w[1]
                 best_model$means[2,2]<-w[2]
                 best_model$means[3,2]<-w[3]
                 best_model$means[4,1]<-w[1]
                 best_model$means[4,2]<-w[4]
                 best_model$covmat[[1]]<-diag(2)
                 best_model$covmat[[2]]<-diag(2)
                 best_model$covmat[[3]]<-diag(2)
                 best_model$covmat[[4]]<-diag(2)
                 best_model$a1=w[5]
                 best_model$a2=w[6]},
           
           mod3={best_model$means[2,1]<-w[1]
                 best_model$means[3,1]<-w[2]
                 best_model$means[3,2]<-w[3]
                 best_model$means[4,1]<-w[4]
                 best_model$means[4,2]<-w[3]
                 best_model$covmat[[1]]<-diag(2)
                 best_model$covmat[[2]]<-diag(2)
                 best_model$covmat[[3]]<-diag(2)
                 best_model$covmat[[4]]<-diag(2)
                 best_model$a1=w[5]
                 best_model$a2=w[6]},
           
           mod4={best_model$means[2,1]<-w[1]
                 best_model$means[3,2]<-w[2]
                 best_model$means[4,1]<-w[1]
                 best_model$means[4,2]<-w[2]
                 vals<-c(1,w[3],w[3],1)
                 best_model$covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$a1=w[4]
                 best_model$a2=w[5]},
          
           
           mod5={best_model$means[2,1]<-w[1]
                 best_model$means[2,2]<-w[2]
                 best_model$means[3,2]<-w[3]
                 best_model$means[4,1]<-w[1]
                 best_model$means[4,2]<-w[4]
                 vals<-c(1,w[5],w[5],1)
                 best_model$covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$a1=w[6]
                 best_model$a2=w[7]},
           
           mod6={best_model$means[2,1]<-w[1]
                 best_model$means[2,2]<-w[2]
                 best_model$means[3,1]<-w[3]
                 best_model$means[3,2]<-w[4]
                 best_model$means[4,1]<-w[5]
                 best_model$means[4,2]<-w[6]
                 best_model$covmat[[1]]<-diag(2)
                 best_model$covmat[[2]]<-diag(2)
                 best_model$covmat[[3]]<-diag(2)
                 best_model$covmat[[4]]<-diag(2)
                 best_model$a1=w[7]
                 best_model$a2=w[8]},
           
           mod7={best_model$means[2,1]<-w[1]
                 best_model$means[3,1]<-w[2]
                 best_model$means[3,2]<-w[3]
                 best_model$means[4,1]<-w[4]
                 best_model$means[4,2]<-w[3]
                 vals<-c(1,w[5],w[5],1)
                 best_model$covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
                 best_model$a1=w[6]
                 best_model$a2=w[7]},
           
           mod8={best_model$means[2,1]<-w[1]
                 best_model$means[3,2]<-w[2]
                 best_model$means[4,1]<-w[1]
                 best_model$means[4,2]<-w[2]
                 best_model$covmat[[1]]<-matrix(c(1,w[3],w[3],1),2,2,byrow=TRUE)
                 best_model$covmat[[2]]<-matrix(c(1,w[4],w[4],1),2,2,byrow=TRUE)
                 best_model$covmat[[3]]<-matrix(c(1,w[5],w[5],1),2,2,byrow=TRUE)
                 best_model$covmat[[4]]<-matrix(c(1,w[6],w[6],1),2,2,byrow=TRUE)
                 best_model$a1=w[7]
                 best_model$a2=w[8]},
           
           mod9={best_model$means[2,1]<-w[1]
                  best_model$means[2,2]<-w[2]
                  best_model$means[3,2]<-w[3]
                  best_model$means[4,1]<-w[1]
                  best_model$means[4,2]<-w[4]
                  best_model$covmat[[1]]<-matrix(c(1,w[5],w[5],1),2,2,byrow=TRUE)
                  best_model$covmat[[2]]<-matrix(c(1,w[6],w[6],1),2,2,byrow=TRUE)
                  best_model$covmat[[3]]<-matrix(c(1,w[7],w[7],1),2,2,byrow=TRUE)
                  best_model$covmat[[4]]<-matrix(c(1,w[8],w[8],1),2,2,byrow=TRUE)
                  best_model$a1=w[9]
                  best_model$a2=w[10]},
           
           mod10={best_model$means[2,1]<-w[1]
                  best_model$means[2,2]<-w[2]
                  best_model$means[3,1]<-w[3]
                  best_model$means[3,2]<-w[4]
                  best_model$means[4,1]<-w[5]
                  best_model$means[4,2]<-w[6]
                  vals<-c(1,w[7],w[7],1)
                  best_model$covmat[[1]]<-matrix(vals,2,2,byrow=TRUE)
                  best_model$covmat[[2]]<-matrix(vals,2,2,byrow=TRUE)
                  best_model$covmat[[3]]<-matrix(vals,2,2,byrow=TRUE)
                  best_model$covmat[[4]]<-matrix(vals,2,2,byrow=TRUE)
                  best_model$a1=w[8]
                  best_model$a2=w[9]},
           
           mod11={best_model$means[2,1]<-w[1]
                  best_model$means[3,1]<-w[2]
                  best_model$means[3,2]<-w[3]
                  best_model$means[4,1]<-w[4]
                  best_model$means[4,2]<-w[3]
                  best_model$covmat[[1]]<-matrix(c(1,w[5],w[5],1),2,2,byrow=TRUE)
                  best_model$covmat[[2]]<-matrix(c(1,w[6],w[6],1),2,2,byrow=TRUE)
                  best_model$covmat[[3]]<-matrix(c(1,w[7],w[7],1),2,2,byrow=TRUE)
                  best_model$covmat[[4]]<-matrix(c(1,w[8],w[8],1),2,2,byrow=TRUE)
                  best_model$a1=w[9]
                  best_model$a2=w[10]},
           
           mod12={best_model$means[2,1]<-w[1]
                  best_model$means[2,2]<-w[2]
                  best_model$means[3,1]<-w[3]
                  best_model$means[3,2]<-w[4]
                  best_model$means[4,1]<-w[5]
                  best_model$means[4,2]<-w[6]
                  best_model$covmat[[1]]<-matrix(c(1,w[7],w[7],1),2,2,byrow=TRUE)
                  best_model$covmat[[2]]<-matrix(c(1,w[8],w[8],1),2,2,byrow=TRUE)
                  best_model$covmat[[3]]<-matrix(c(1,w[9],w[9],1),2,2,byrow=TRUE)
                  best_model$covmat[[4]]<-matrix(c(1,w[10],w[10],1),2,2,byrow=TRUE)
                  best_model$a1=w[11]
                  best_model$a2=w[12]})
    
    # save convergence info for best-fitting model
    best_model$convergence <- optim_list[[model_num]]$convergence
    best_model$message <- optim_list[[model_num]]$message
    
    # get observed and predicted values
    best_model$predicted <- as.vector(matrix_predict(best_model$means, best_model$covmat, diag(2), matrix(c(best_model$a1,best_model$a2),2,1)))
    best_model$observed <- as.vector(pmatrix(data))
    
    # return object of class grt_hm_fit
    results <- list(table=o_aic, best_model=best_model)
    class(results) <- "grt_hm_fit"
    return (results)
}
