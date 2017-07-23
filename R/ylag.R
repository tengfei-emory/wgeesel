ylag<-function(id,y,lag,na=FALSE){
  res=lapply(unique(id),function(x){
     idx=which(id==x)
     y_i=y[idx]
     c(rep(NA,lag),y_i[1:(length(y_i)-lag)])
     })
  res=unlist(res)
  if (lag>1 & na==FALSE){
    res= ifelse(is.na(res),0,res)
  }
  return(res)
}


