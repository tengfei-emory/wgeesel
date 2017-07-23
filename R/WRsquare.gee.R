WRsquare.gee<-function(object,weight_mean){
  #fit=wgee(model=model,data=data,id=id,family=family,corstr=corstr,scale=scale,mismodel=mismodel)
  mu_fit=object$mu_fit
  weight=object$weight
  y=object$y
  if(weight_mean==T) Y_bar=sum(weight*y,na.rm=T)/length(y)
  else Y_bar=mean(y,na.rm=T)
  WRsquare=sum(weight*(y-mu_fit)^2,na.rm=T)/sum(weight*(y-Y_bar)^2,na.rm=T)
  return(list(WRsquare=WRsquare))
}
