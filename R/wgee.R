wgee<-function(model,data,id,family,corstr,scale=NULL,mismodel=NULL){
  call <- match.call()
  m=model.frame(model,data,na.action="na.pass")
  y=model.response(m,"numeric")
  x=model.matrix(model,m)
  if(!is.null(mismodel)){
    rm=model.frame(mismodel,data,na.action="na.pass")
    r=model.response(rm,"numeric")
    z=model.matrix(mismodel,rm)
    adjusted_idx=lapply(1:length(unique(id)),function(x){
      idx=which(id==unique(id)[x])
      mis_idx=which(is.na(y[idx])==T)
      if (length(mis_idx)>0) a=idx[2:min(mis_idx)]
      else a=idx[2:length(idx)]
      a
    })
    adjusted_idx=unlist(adjusted_idx)
    
    mis_fit=glm(mismodel,data=data[adjusted_idx,],family=binomial())
    
    adjusted_w=lapply(1:length(unique(id)),function(x){
      idx=which(id==unique(id)[x])
      predict_d=data[idx,]
      predict_w=predict(mis_fit,newdata=data[idx,],type="response")
      predict_w[1]=1
      predict_w=unlist(lapply(1:length(idx),function(m){prod(predict_w[1:m])}))
      res=r[idx]/predict_w
      res[which(is.na(res))]=0
      return(res)
    })
    weight=unlist(adjusted_w)
  
  fit=newton_raphson(id,x,y,weight=weight,scale=scale,corstr=corstr,family=family)

  beta_est=fit$beta_new
  rownames(beta_est)=colnames(x)

  scale=fit$phi
  rho=fit$rho
  res_list=lapply(1:length(unique(id)),function(m){
    idx=which(id==unique(id)[m])
    us_matrix(x[idx,],y[idx],weight[idx],beta_est,rho,scale,corstr,family,z[idx,],r[idx],coef(mis_fit))
  })
  U_i=lapply(res_list,function(x){x[[1]]})
  logit_S_i=lapply(res_list,function(x){x[[2]]})
  US_i=lapply(res_list,function(x){x[[3]]})
  SS_i=lapply(res_list,function(x){x[[4]]})
  D_i=lapply(res_list,function(x){x[[5]]})
  W_i=lapply(res_list,function(x){x[[6]]})
  V_i=lapply(res_list,function(x){x[[7]]})
  sum_US_SS_i=Reduce("+",US_i)%*%solve(Reduce("+",SS_i))
  variance=V_w_est(id,U_i,logit_S_i,sum_US_SS_i,D_i,W_i,V_i)
  mu_fit=exp(x%*%beta_est)/(1+exp(x%*%beta_est))
  final_res=list(beta=beta_est,var=variance,mu_fit=mu_fit,scale=scale,rho=fit$rho,weight=weight,model=model,x=x,y=y,mis_fit=mis_fit,call=call,id=id,data=data,family=family,corstr=corstr)
  }
  if(is.null(mismodel)){
   
   data$id=id
   fit=geeglm(model,data=data,id=id,family=family,corstr = corstr,scale.fix = !is.null(scale))
   final_res=list(beta=fit$geese$beta,var=fit$geese$vbeta,mu_fit=fit$fitted.values,scale=fit$geese$gamma,rho=fit$geese$alpha,weight=rep(1,nrow(data)),model=model,x=x,y=y,mis_fit=NA,call=call,id=id,data=data,family=family,corstr=corstr)
   
  }
  class(final_res)=c("wgee") 
  return(final_res)
}

print.wgee <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(x$beta))
  cat("\n Scale Parameter: ", signif(x$scale, digits=4), "\n")
  cat("\n Estimated Correlation Parameter: ", signif(x$rho, digits=4), "\n")
}

summary.wgee<- function(object, ...)  {
  Coefs <- matrix(NA,nrow=length(object$beta),ncol=4)
  Coefs[,1] <- c(object$beta)
  Coefs[,2] <- sqrt(diag(object$var))
  Coefs[,3] <- Coefs[,1]/Coefs[,2]
  Coefs[,4] <- round(2*pnorm(abs(Coefs[,3]), lower.tail=F), digits=8)
  colnames(Coefs) <- c("Estimates","Robust SE", "z value", "Pr(>|z|)")
  coefnames<-rownames(object$beta)
  summ <- list(beta = Coefs[,1],se.robust = Coefs[,2], wald.test = Coefs[,3], p = Coefs[,4],
               corr = object$rho, phi = object$scale, call=object$call,coefnames=coefnames)
  class(summ) <- 'summary.wgee'
  return(summ)
}

print.summary.wgee<- function(x,digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  Coefs <- matrix(0,nrow=length(x$beta),ncol=4)
  rownames(Coefs) <- c(x$coefnames)
  colnames(Coefs) <- c("Estimates","Robust SE", "z value", "Pr(>|z|)")
  Coefs[,1] <- x$beta
  Coefs[,2] <- x$se.robust
  Coefs[,3] <- x$wald.test
  Coefs[,4] <- x$p
  printCoefmat(as.matrix(Coefs), digits = digits) ## Thanks, Achim
  cat("\n Estimated Correlation Parameter: ", signif(x$corr, digits=4), "\n")
  cat("Estimated Scale Parameter: ", signif(x$phi, digits=4), "\n")
}

