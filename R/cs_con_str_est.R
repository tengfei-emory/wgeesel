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
  e_i=(y_i-mu_i)/sqrt(cons_cov)
  mat<-diag(w_i)
  big<-e_i%*%t(e_i)  #what i change
  big[is.na(big)]<-0  #what i change
  std_resi=big%*%mat   #what i change
  #std_resi=(w_i*e_i)%*%t(w_i*e_i)
  #std_resi[is.na(std_resi)]<-0
  res=0
  if(corstr=="ar1"){
   res=sum(diag(std_resi[1:(nrow(std_resi)-1),2:ncol(std_resi)]),na.rm=T)
  }
  else if(corstr=="unstructured"){
    std_resi[lower.tri(std_resi,diag=F)]=std_resi[upper.tri(std_resi,diag=F)]
    res=std_resi
    #res=big
  }
  else{
    res=sum(std_resi[upper.tri(std_resi,diag=F)],na.rm=T)
    }
  res
}
