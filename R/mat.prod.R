mat.prod <-
function(A,B){
    if (is.matrix(A)==FALSE) {A=t(A)
    }else{A=A}
   matrices <- matrix(0, ncol=ncol(A), nrow=nrow(A))
    for (i in 1:ncol(A)) {
      matrices[,i] <- A[,i]*B
    }
   return(matrices)
  }
