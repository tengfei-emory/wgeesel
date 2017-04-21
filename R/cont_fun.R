cont_fun <- function(beta,rho,phi,x,id, corstr){ #N is the sample size, nsubj is cluster size;
  response=c()
  N=length(unique(id))
  for(i in 1:N){
    mu_i <- x[which(id==i),]%*%beta 
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
    
    S_i=diag(rep(1,length(mu_i)))
    m=phi*(sqrt(S_i)%*%R_i%*%sqrt(S_i))
    response <- c(response,mvrnorm(n = 1, mu=mu_i, Sigma=m))
  }
  data=cbind(id,x,response)
  data=data.frame(data)
  
  return(data)
}