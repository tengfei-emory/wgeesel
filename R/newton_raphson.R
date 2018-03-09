newton_raphson <-function(id,x,y,weight,scale,corstr,family,maxit,tol){
  id_uni=unique(id)
  time<-length(id)/length(id_uni)
  d_obs=na.omit(data.frame(cbind(y,x,id)))
  #d_obs=data.frame(cbind(y,x,id))
  beta_s=coef(glm(y~x-1,family = family))
  n_i=table(d_obs$id)
  N=nrow(d_obs)
  if(is.null(weight)){
    p1=length(beta_s)
    p2=length(beta_s)
    weight=rep(1,N)
    }
  else {
    p1=0
    p2=0
  }
  iter<-0
  repeat{
    if(iter>maxit) {
      warning("Please ncrease maximum iteration number, the algorithm is not converged")
      break
    }
    phi_vector=lapply(id_uni,function(m){
      idx=which(id==m)
      phi_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
    )
    if(is.null(scale)) {
      phi=sum(unlist(phi_vector))/(length(id)-p1) #change denorminator
    }
    else {phi=scale}
    
    if(corstr=="independence") 
      {rho=0}
    else if(corstr=="exchangeable"){
      rho_vector=lapply(id_uni,function(m){
      idx=which(id==m)
      cs_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
      )
      #N_star=0.5*sum(n_i*(n_i-1))
      N_star=0.5*length(id_uni)*(time-1)*time #I change here
      rho=sum(unlist(rho_vector))/(phi*(N_star-p2))
    }
    else if(corstr=="ar1") {
      rho_vector=lapply(id_uni,function(m){
        idx=which(id==m)
        cs_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
      )
      #N_star=sum(n_i-1)
      N_star<-length(id_uni)*(time-1) #here I change 
      rho=sum(unlist(rho_vector))/(phi*(N_star-p2))
    }
    else if(corstr=="unstructured") {
      rho_vector=lapply(id_uni,function(m){
        idx=which(id==m)
        cs_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
      )

      rho=Reduce("+",rho_vector)/(phi*(length(id_uni)-p2))
      rho=rho-diag(diag(rho),nrow(rho))+diag(1,nrow(rho))
    }
    
    est=lapply(id_uni,function(m){
      idx=which(id==m)
      cluster_est(x[idx,],y[idx],weight[idx],beta_s,rho,phi,corstr,family)
    })
  
    mat1=lapply(est,function(x){x[[1]]})
    mat2=lapply(est,function(x){x[[2]]}) 

    beta_new=beta_s+solve(Reduce('+',mat2))%*%Reduce('+',mat1)
    dif=sum(abs(beta_new-beta_s))
    if(dif<tol) {
      beta_s=beta_new
      iter<-iter+1 
      break
    }
    else {
      beta_s=beta_new
      iter<-iter+1 
    }
  }
  list(beta_new=beta_new,rho=rho,phi=phi)
}
