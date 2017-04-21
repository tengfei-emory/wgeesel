V_w_est <-
function(id,U,logit_S,sumUS,D,W,V){
  mat=lapply(1:length(unique(id)),function(m_v){
    E_i_est(m_v,U,logit_S,sumUS,D,W,V)
  })
 x_mat_1=Reduce("+",lapply(mat,function(x){x[[1]]}))
 x_mat_2=Reduce("+",lapply(mat,function(x){x[[2]]}))
 V_w=solve(x_mat_2)%*%x_mat_1%*%solve(x_mat_2)
 return(V_w)
}
