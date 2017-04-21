binary_fun <- function(beta,rho,x,id, corstr){ #N is the sample size, nsubj is cluster size;
  response=c()
  N=length(unique(id))
  for(i in 1:N){
    mu_i <- 1/(1+exp(-(x[which(id==i),]%*%beta)))
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
    response <- c(response, rmvbin(1, margprob = as.vector(mu_i), bincorr = R_i))
  }
  data=cbind(id,x,response)
  data=data.frame(data)
  return(data)
}