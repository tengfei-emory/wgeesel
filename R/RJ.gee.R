RJ.gee <-
function(object){
  model=object$model
  data=object$data
  id=object$id
  corstr=object$corstr
  family=object$family
  
  ####beta is model estimates###
  if (length(id) != nrow(data)){
    stop("variable lengths differ (found for '(id)')")}
  m.temp <- model.frame(model, data, na.action='na.pass')
  cluster <- cluster.size(id)$n #number of observations for each subject;
  size <- cluster.size(id)$m # sample size;
  subject <- rep(1:size, cluster) #1:size; create subject variable 
  
  data.temp <- cbind(subject=subject,m.temp) #combine the design matrix and subject varible
  ### get rid of the missing data;
  data <- na.omit(data.temp)
  subject=data$subject
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
  step11<-matrix(0, nrow=len, ncol=len)
  for (i in 1:size){
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(subject)[i]))
    var_i<-var[1:cluster[i],1:cluster[i]]
    nsubj_i <- cluster[i]
    switch(
      family,
      gaussian={ xx<-t(covariate)%*%solve(var_i)%*%covariate
      },
      binomial={ u<-1/as.vector(1+exp(-covariate%*%beta))
      D<-covariate*u
      diag <-apply(as.matrix(u), 1, function(x) {sqrt(x*(1-x))})
      A<- solve(diag(diag,nsubj_i))
      xx<-t(D)%*%A%*%solve(var_i)%*%A%*%D
      },
      poisson={ B<-matrix(0,nrow=nsubj_i,ncol=nsubj_i)
      diag(B)<-1/sqrt(exp(covariate%*%beta))
      #D<-cbind(exp(covariate%*%beta), m1*exp(covariate%*%beta), m2*exp(covariate%*%beta))
      D <- mat.prod(covariate, exp(covariate%*%beta))
      xx<-t(D)%*%B%*%solve(var_i)%*%B%*%D
      }
    )
    step11<-step11+xx
  }     
  
  step12<-matrix(0,nrow=len,ncol=len)
  for (i in 1:size){
    nsubj_i <- cluster[i]
    y_i<-as.matrix(data[subject==i,]$response)
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(subject)[i]))
    var_i<-var[1:cluster[i],1:cluster[i]]
    switch(
      family,
      gaussian={ xy<-t(covariate)%*%solve(var_i)%*%(y_i-covariate%*%beta)
      },
      binomial={  u<-1/as.vector(1+exp(-covariate%*%beta))
      D<-covariate*u
      diag <-apply(as.matrix(u), 1, function(x) {sqrt(x*(1-x))})
      A<- solve(diag(diag,nsubj_i))
      xy<-t(D)%*%A%*%solve(var_i)%*%A%*%(y_i -1/(1+exp(-covariate%*%beta)))
      },
      poisson={ B<-matrix(0,nrow=nsubj_i,ncol=nsubj_i)
      diag(B)<-1/sqrt(exp(covariate%*%beta))
      D <- mat.prod(covariate, exp(covariate%*%beta))
      #D<-cbind(exp(covariate%*%beta), m1*exp(covariate%*%beta), m2*exp(covariate%*%beta))
      xy<-t(D)%*%B%*%solve(var_i)%*%B%*%(y_i -exp(covariate%*%beta))
      }
    )
    step12<-step12+xy%*%t(xy)
  }
  
  Q<-solve(step11)%*%(step12)
  c1<-sum(diag(Q))/3
  c2<-sum(diag(Q%*%Q))/3
  rj<-sqrt((1-c1)^2+(1-c2)^2)
  return(list(RJC=rj))
}
