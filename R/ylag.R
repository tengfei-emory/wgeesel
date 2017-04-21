ylag<-function(id,y,lag){
  res=lapply(unique(id),function(x){
     idx=which(id==x)
     y_i=y[idx]
     c(rep(NA,lag),y_i[1:(length(y_i)-lag)])
     })
  res=unlist(res)
  return(res)
}


