cs_con_str_est <-
function(x_i,y_i,w_i,beta,corstr,family){
  if(family=="binomial"){
    mu_i=exp(x_i%*%beta)/(1+exp(x_i%*%beta))
    cons_cov=c(mu_i*(1-mu_i))
  }
  else if (family=="gaussian") {
    mu_i=x_i%*%beta
    cons_cov=rep(1,length(mu_i))
  }
  else {
    mu_i=exp(x_i%*%beta)
    cons_cov=c(mu_i)
  } 
  e_i=w_i*(y_i-mu_i)/sqrt(cons_cov)
  std_resi=e_i%*%t(e_i)
  res=0
  if(corstr=="ar1"){
   res=sum(diag(std_resi[1:(nrow(std_resi)-1),2:ncol(std_resi)]),na.rm=T)
  }
  else{
    res=sum(std_resi[upper.tri(std_resi,diag=F)],na.rm=T)
    }
  res
}
