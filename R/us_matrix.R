us_matrix=function(x_i,y_i,w_i,beta,rho,phi,corstr,family,z_i,r_i,alpha){
  mat_est=list()
  non_missing_length=length(na.omit(y_i))
  if(family=="binomial"){
    mu_i=exp(x_i%*%beta)/(1+exp(x_i%*%beta))
    cons_cov=c(mu_i*(1-mu_i)) #n_i*1 matrix#
    D_i=x_i*cons_cov #n_i*p#
  }
  else if (family=="gaussian"){
    mu_i=x_i%*%beta
    cons_cov=rep(1,length(mu_i)) #n_i*1 matrix#
    D_i=x_i #n_i*p#
  }
  else if (family=="poisson"){
    mu_i=exp(x_i%*%beta)
    cons_cov=c(mu_i) #n_i*1 matrix#
    D_i=x_i*cons_cov #n_i*p#
  }
  S_i=diag(cons_cov) #n_i*n_i#
  W_i=diag(w_i)
  if(corstr=="exchangeable"){
    R_i=diag(1,length(mu_i))+matrix(rho,ncol=length(mu_i),nrow=length(mu_i))-diag(rho,length(mu_i))
  }
  if(corstr=="independence"){
    R_i=diag(1,length(mu_i))
  }
  if(corstr=="ar1"){
    H <- abs(outer(1:length(mu_i),1:length(mu_i), "-")) 
    R_i=rho^H
  }
  if(corstr=="unstructured"){
    R_i=rho
  }
  V_i=phi*(sqrt(S_i)%*%R_i%*%sqrt(S_i))
  U_i=(t(D_i)%*%solve(V_i)%*%W_i)[,1:non_missing_length]%*%as.matrix(((y_i-mu_i)[1:non_missing_length]),ncol=1)
  lam_it=exp(z_i%*%alpha)/(1+exp(z_i%*%alpha))
  
  logit_S_i_list=lapply(2:length(r_i),function(x){
    res=0
    if(r_i[x-1]==1)res=(r_i[x-1])*(r_i[x]-lam_it[x])*t(t(z_i[x,]))
    return(res)
  })
  logit_S_i=Reduce("+",logit_S_i_list)
  mat_est[[1]]=U_i
  mat_est[[2]]=logit_S_i
  mat_est[[3]]=U_i%*%t(logit_S_i)
  mat_est[[4]]=logit_S_i%*%t(logit_S_i)
  mat_est[[5]]=D_i
  mat_est[[6]]=W_i
  mat_est[[7]]=V_i
  return(mat_est)
}