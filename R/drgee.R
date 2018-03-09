drgee<-function(model,outcomemodel,data,id,family,corstr,scale=NULL,mismodel=NULL,nameTRT,maxit=200, tol=0.001){
    call <- match.call()
    m=model.frame(model,data,na.action="na.pass")
    y=model.response(m,"numeric")
    x=model.matrix(model,m)
    if(!is.null(mismodel)){
      rm=model.frame(mismodel,data,na.action="na.pass")
      r=model.response(rm,"numeric")
      z=model.matrix(mismodel,rm)
      adjusted_idx=lapply(1:length(unique(id)),function(x){
        idx=which(id==unique(id)[x])
        mis_idx=which(is.na(y[idx])==T)
        if (length(mis_idx)>0) a=idx[2:min(mis_idx)]
        else a=idx[2:length(idx)]
        a
      })
      adjusted_idx=unlist(adjusted_idx)
      
      mis_fit=glm(mismodel,data=data[adjusted_idx,],family=binomial())
      
      adjusted_w=lapply(1:length(unique(id)),function(x){
        idx=which(id==unique(id)[x])
        predict_d=data[idx,]
        predict_w=predict(mis_fit,newdata=data[idx,],type="response")
        predict_w[1]=1
        predict_w=unlist(lapply(1:length(idx),function(m){prod(predict_w[1:m])}))
        res=r[idx]/predict_w
        res[which(is.na(res))]=0
        return(res)
      })
      weight=unlist(adjusted_w)
    }
    
    weight=rep(1,nrow(data))
      fit=geeDREstimation(formula=model,
                          id=id , data = data,
                          family = family, corstr = corstr, weights=ifelse(weight==0,0,1/weight),
                          model.augmentation.trt=outcomemodel,
                          model.augmentation.ctrl=outcomemodel,typeweights = "VW",
                          scale.fix=!is.null(scale),init.phi=ifelse(!is.null(scale),scale,1),nameTRT=nameTRT,maxit=maxit,tol=tol)
     return(fit)
}


