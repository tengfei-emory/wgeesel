newton_raphson <-function(id,x,y,weight,scale,corstr,family){
  id_uni=unique(id)
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
  else{
    p1=0
    p2=0
  }
  repeat{
    phi_vector=lapply(id_uni,function(m){
      idx=which(id==m)
      phi_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
    )
    if(is.null(scale)) {
      phi=sum(unlist(phi_vector))/(N-p1)
    }
    else {phi=scale}
    
    if(corstr=="independence") {rho=0}
    else if(corstr=="exchangeable"){
      rho_vector=lapply(id_uni,function(m){
      idx=which(id==m)
      cs_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
      )
      N_star=0.5*sum(n_i*(n_i-1))
      rho=sum(unlist(rho_vector))/(phi*(N_star-p2))
    }
    else if(corstr=="ar1") {
      rho_vector=lapply(id_uni,function(m){
        idx=which(id==m)
        cs_con_str_est(x[idx,],y[idx],weight[idx],beta_s,corstr,family)}
      )
      N_star=sum(n_i-1)
      rho=sum(unlist(rho_vector))/(phi*(N_star-p2))
    }
    est=lapply(id_uni,function(m){
      idx=which(id==m)
      cluster_est(x[idx,],y[idx],weight[idx],beta_s,rho,phi,corstr,family)
    })
  
    mat1=lapply(est,function(x){x[[1]]})
    mat2=lapply(est,function(x){x[[2]]}) 
    beta_new=beta_s+solve(Reduce('+',mat2))%*%Reduce('+',mat1)
    dif=sum(abs(beta_new-beta_s))
    if(dif<0.001) {
      break
    }
    else {
      beta_s=beta_new
     }
   }
  list(beta_new=beta_new,rho=rho,phi=phi)
}
