RJ2.gee <-
function(object){
  
  model=object$model
  mismodel=object$mismodel
  data=object$data
  id=object$id
  corstr=object$corstr
  family=object$family
  
  if (sum(is.na(id)) !=0 ){
    stop("ID has missing")}
  
  m.temp <- model.frame(model, data, na.action='na.pass')
  cluster <- cluster.size(id)$n
  ###check if the model has missing
  if (sum(is.na(m.temp)) !=0 | max(cluster) != min(cluster)){
    stop("The data has missing, cluster size is not equal")}
  
  if (length(id) != nrow(m.temp)){
    stop("variable lengths differ (found for '(id)')") }
  
  #number of observations for each subject;
  size <- cluster.size(id)$m # sample size;
  subject <- rep(1:size, cluster) #1:size; create subject variable 
  
  data <- cbind(subject,m.temp) #combine the design matrix and subject varible
  
  ### Get the design matrix after removing the missing data
  m <- model.frame(model, data, na.action='na.pass')
  data$response <- model.response(m, "numeric")
  cluster <- cluster.size(subject)$n #number of observations for each subject;
  nsubj <- max(cluster)
  size <- cluster.size(subject)$m # sample size;
  mat <- as.data.frame(model.matrix(model, m)) #covariate matrix with intercept
  mat$subj <- rep(unique(subject), cluster)
  
  switch(family,
         gaussian={model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=gaussian)
         },
         binomial={model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=gaussian)
         },
         poisson={model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=gaussian)
         },
         stop("Warnings: Invalid type of outcomes!")
  )
  beta <- as.matrix(model.R$coef)  
  scale <- as.numeric(summary(model.R)$geese$scale[1])
  var<- switch(corstr,
               independence=scale*cormax.ind(nsubj),
               exchangeable=scale*cormax.exch(nsubj, model.R$geese$alpha), ###scale*cormax.ind
               ar1=scale*cormax.ar1(nsubj, model.R$geese$alpha)
  )
  len <- length(beta)
  step01<-matrix(0, nrow=len, ncol=len)
  for (i in 1:size){
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(subject)[i]))
    var_i<-var[1:cluster[i],1:cluster[i]]
    switch(
      family,
      gaussian={ xx<-t(covariate)%*%solve(var_i)%*%covariate
      },
      binomial={ u<-1/as.vector(1+exp(-covariate%*%beta))
      D<-covariate*u
      diag <-apply(as.matrix(u), 1, function(x) {sqrt(x*(1-x))})
      A<- solve(diag(diag,nsubj))
      xx<-t(D)%*%A%*%solve(var_i)%*%A%*%D
      },
      poisson={ 
        B<-matrix(0,nrow=nsubj,ncol=nsubj)
        diag(B)<-1/sqrt(exp(covariate%*%beta))
        D <- mat.prod(covariate, exp(covariate%*%beta))
        xx<-t(D)%*%B%*%solve(var_i)%*%B%*%D
      }
    )
    step01<-step01+xx
  }     
  step<-matrix(0, nrow=nsubj, ncol=nsubj)
  for (i in 1:size){
    y_i<-as.matrix(data[subject==i,]$response)
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(subject)[i]))
    var_i<-var[1:cluster[i],1:cluster[i]]
    switch(
      family,
      gaussian={ resid<-solve(cormax.ind(nsubj)-covariate%*%solve(step01)%*%t(covariate)
                              %*%solve(var_i))%*%(y_i-covariate%*%beta)
      },
      binomial={  u<-1/as.vector(1+exp(-covariate%*%beta))
      D<-covariate*u
      diag<-apply(as.matrix(u), 1, function(x) {sqrt(x*(1-x))})
      A<-solve(diag(diag,nsubj))
      resid<-A%*%solve(cormax.ind(nsubj)-D%*%solve(step01)%*%t(D)%*%A%*%solve(var_i)%*%A)%*%(y_i-1/(1+exp(-covariate%*%beta)))
      },
      poisson={ B <- matrix(0,nrow=nsubj,ncol=nsubj)
      diag(B)<-1/sqrt(exp(covariate%*%beta))
      D <- mat.prod(covariate, exp(covariate%*%beta))
      resid<-B%*%solve(cormax.ind(nsubj)-D%*%solve(step01)%*%t(D)%*%B%*%solve(var_i)%*%B)%*%(y_i-exp(covariate%*%beta))
      }
    )
    step<-step+resid%*%t(resid)
  }
  #unstr <- matrix(0,nrow=len,ncol=len) #### same for each subject??? cov(y_i)
  unstr<-step/size
  diag(unstr)<-rep(1, nsubj)
  step11<-matrix(0,nrow=len,ncol=len)
  step12<-matrix(0,nrow=len,ncol=len)
  for (i in 1:size){
    y_i <- as.matrix(data[subject==i,]$response)
    covariate <- as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(subject)[i]))
    var_i <- var[1:cluster[i],1:cluster[i]]
    switch(
      family,
      gaussian={  xy<-t(covariate)%*%solve(var_i)%*%unstr%*%solve(var_i)%*%covariate
      xx<-t(covariate)%*%solve(var_i)%*%covariate
      },
      binomial={  u<-1/as.vector(1+exp(-covariate%*%beta))
      D<-covariate*u
      diag<-apply(as.matrix(u), 1, function(x) {sqrt(x*(1-x))})
      A<-solve(diag(diag,nsubj))
      xy<-t(D)%*%A%*%solve(var_i)%*%unstr%*%solve(var_i)%*%A%*%D
      xx<-t(D)%*%A%*%solve(var_i)%*%A%*%D
      },
      poisson={ B<-matrix(0,nrow=nsubj,ncol=nsubj)
      diag(B)<-sqrt(exp(covariate%*%beta))
      D <- mat.prod(covariate, exp(covariate%*%beta))
      xy<-t(D)%*%solve(B)%*%solve(var_i)%*%unstr%*%solve(var_i)%*%solve(B)%*%D
      xx<-t(D)%*%solve(B)%*%solve(var_i)%*%solve(B)%*%D
      }
    )
    step11<-step11+xx
    step12<-step12+xy
  }
  
  Q<-solve(step11)%*%(step12)
  c1<-sum(diag(Q))/3
  c2<-sum(diag(Q%*%Q))/3
  rj<-sqrt((1-c1)^2+(1-c2)^2)
  return(list(RJC=rj))
}
