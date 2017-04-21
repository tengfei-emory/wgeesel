wgee<-function(model,data,id,family,corstr,scale=NULL,mismodel=NULL){
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
  }
  if(is.list(mismodel)){weight=rep(1,length(y))}
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
  se=sqrt(diag(variance))
  p_value=pchisq((beta_est/se)^2,1,lower.tail = F)
  print(data.frame(est=beta_est,se=se,p_value=round(p_value,4)))
  final_res=list(beta=beta_est,var=var,mu_fit=mu_fit,scale=scale,rho=fit$rho,weight=weight,model=model,x=x,y=y,mis_fit=mis_fit)
  return(final_res)
}