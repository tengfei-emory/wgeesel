cluster_est <-
function(x_i,y_i,w_i,beta,rho,phi,corstr,family){
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
    cons_cov=c(exp(x_i%*%beta)) #n_i*1 matrix#
    D_i=x_i*cons_cov #n_i*p#
  }
  S_i=diag(cons_cov) #n_i*n_i#
  W_i=diag(w_i)
  if(corstr=="exchangeable"){
    R_i=diag(1,length(mu_i))+matrix(rho,ncol=length(mu_i),nrow=length(mu_i))-diag(rho,length(mu_i))
  }
  else if(corstr=="independence"){
    R_i=diag(1,length(mu_i))
  }
  else if(corstr=="ar1"){
    H <- abs(outer(1:length(mu_i),1:length(mu_i), "-")) 
    R_i=rho^H
  }
  V_i=phi*(sqrt(S_i)%*%R_i%*%sqrt(S_i))
  mat_est[[1]]=(t(D_i)%*%solve(V_i)%*%W_i)[,1:non_missing_length]%*%as.matrix((y_i-mu_i)[1:non_missing_length],ncol=1)
  mat_est[[2]]=t(D_i)%*%solve(V_i)%*%D_i
  return(mat_est)
}
