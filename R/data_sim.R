data_sim <- function(id, rho,phi, x, beta, x_mis, para, corstr, family,lag_level){
 
  switch(family,
         "gaussian"={data <- cont_fun(beta,rho,phi,x,id,corstr)
         },
         "binary"={data <- binary_fun(beta,rho,x,id,corstr)
         },
         "poisson"={data <- count_fun(beta,rho,x,id,corstr)
         },
         stop("Warnings: Invalid type!")
  )
  complete=F
  x_final=list()
  prob_miss=0
  if(complete==F){
    data$ind <- 1
    #data$prey_mis <- 0
    data$y_first <- NULL;
    N=length(unique(id))
    for (i in 1:N){
      nsubj=max(table(id))
      for (k in 2:nsubj){
      lagy=NULL
      if(lag_level!=0){
      for(lag_i in 1:lag_level){
           ylagname=paste("ylag",lag_i,sep="")
           y_i=data$response[which(id==i)]
           y_new=c(rep(NA,lag_i),y_i[1:(length(y_i)-lag_i)])
           lagy=cbind(lagy,y_new)
           colnames(lagy)[lag_i]=ylagname
      }
      }
       
        x_mis_i=cbind(x_mis[which(id==i),],lagy=lagy)
        x_final[[i]]=x_mis_i
         
        p <- 1/(1+exp(sum(-x_mis_i[k,]*para,na.rm=T)))
        ###depends on the first and second observation;
        ind_p <- rbinom(1,1,p)
        data[data$id==i,]$ind[k] <- ind_p 
        if (data[data$id==i,]$ind[k-1]==0){
          data[data$id==i,]$ind[(k):nsubj] <- rep(0,(nsubj-(k)+1))
        }
      }  
      data$y_first[(i*nsubj-(nsubj-1)):(i*nsubj)]=rep(data[data$id==i,]$response[1],nsubj)
    }
    x_final=do.call("rbind",x_final)
    response_mis=ifelse(data$ind==0,NA,data$response)
    data_final=cbind(data[,c("id","response","ind")],response_mis,x,x_final)
    data<-as.data.frame(data_final)
    prob_miss<-length(which(data$ind==0))/length(data$response)
  }
  return(list(data=data, prob_miss=prob_miss))
}
