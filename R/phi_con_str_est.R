phi_con_str_est <-
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
  adj_factor=w_i
  e_i=adj_factor*(y_i-mu_i)/sqrt(cons_cov)
  std_resi=e_i%*%t(e_i)
  sum(c(diag(std_resi)),na.rm=T)
}
