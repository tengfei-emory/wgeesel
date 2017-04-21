E_i_est <-
function(id_e_i,U_i,logit_S_i,sumUS,D_i,W_i,V_i){
  mat_est=list()                  
  mat1=U_i[[id_e_i]]-sumUS%*%logit_S_i[[id_e_i]]
  mat_est[[1]]=mat1%*%t(mat1)
  mat_est[[2]]=t(D_i[[id_e_i]])%*%solve(V_i[[id_e_i]])%*%W_i[[id_e_i]]%*%D_i[[id_e_i]]
  return(mat_est)
}
