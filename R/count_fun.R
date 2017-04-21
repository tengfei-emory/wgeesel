count_fun <- function(beta,rho,x,id, corstr){ 
  response=c()    
  N=length(unique(id))
  for(i in 1:N){
    mu_i <- exp(x[which(id==i),]%*%beta)
    if(corstr=="independence"){
      R_i=diag(1,length(mu_i))
    }
    else if (corstr=="exchangeable"){
      R_i=diag(1,length(mu_i))+matrix(rho,ncol=length(mu_i),nrow=length(mu_i))-diag(rho,length(mu_i))
    }
    else {
      H <- abs(outer(1:length(mu_i),1:length(mu_i), "-")) 
      R_i=rho^H
    }
    y <- genPoisNor(n=1, no.pois=length(mu_i), cmat.star=R_i,
                    lamvec=mu_i,no.norm=0,mean.vec=NULL,sd.vec=NULL)
    response <- c(response, y)
  }
  data=cbind(id,x,response)
  data=data.frame(data)
  return(data)
}
