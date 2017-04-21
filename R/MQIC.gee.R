MQIC.gee <-
function(model, data, id,family, corstr) {
  ### data information;
  if (length(id) != nrow(data)){
    stop("variable lengths differ (found for '(id)')")}  
  m1 <- model.frame(model, data, na.action='na.pass')
  m1$id <- id
  data <- na.omit(m1)

  cluster<-cluster.size(data$id)$n #number of observations for each subject;
  size<-cluster.size(data$id)$m # sample size;
  browser()
  data$subject <- rep(1:size, cluster) #1:size;
  
  ### Get the design matrix;
  m <- model.frame(model, data, na.action='na.pass')
  data$response <- model.response(m, "numeric")
  mat <- as.data.frame(model.matrix(model, m))
  mat$subj <- rep(unique(data$subject), cluster)
  
  subject <- data$subject
  # Quasi Likelihood;
  switch(family,
         gaussian={model.ind <- geeglm(model, id =subject, data=data, corstr="independence", family=gaussian)
         ### alternative correlation structure;
         model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=gaussian)
         beta <- as.matrix(as.numeric(model.R$coefficients))
         scale <- as.numeric(summary(model.R)$geese$scale[1])
         mu.R <- model.R$fitted.values
         data$mu.R <- mu.R
         y <- model.R$y
         quasi.R <- -1/2*sum((y-mu.R)^2)     
         },
        
         binomial={model.ind <- geeglm(model, id =subject, data=data, corstr="independence", family=binomial)
         ### alternative correlation structure;
         model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=binomial)
         beta <- as.matrix(as.numeric(model.R$coefficient))
         scale <- as.numeric(summary(model.R)$geese$scale[1])
         mu.R <- model.R$fitted.values
         data$mu.R <- mu.R
         y <- model.R$y
         quasi.R <- sum(y*log(mu.R/(1-mu.R))+log(1-mu.R))
         },
         
         poisson={model.ind <- geeglm(model, id =subject, data=data, corstr="independence", family=poisson)
         ### alternative correlation structure;
         model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=poisson)
         beta <- as.matrix(as.numeric(model.R$coefficients))
         scale <- as.numeric(summary(model.R)$geese$scale[1])
         mu.R <- model.R$fitted.values
         data$mu.R <- mu.R
         y <- model.R$y
         quasi.R <- sum(y*log(mu.R)-mu.R)
         },
         stop("Warnings: Invalid type of outcomes!")
  )
  Vnaive <- model.R$geese$vbeta.naiv
  px <- dim(Vnaive)[1] # number non-redunant columns in design matrix
  
  # calculate Va=\sum(Di%*%Vi(-1)%*%Cov(y)%*%A^{-1}%*%Di);
  va<-matrix(0,nrow=px,ncol=px) #p by p matrix
  for (i in 1:size){
    ncluster <- cluster[i]
    y_i<-as.matrix(data[data$subject==i,]$response)
    mu.R_i <-as.matrix(data[data$subject==i,]$mu.R) 
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$subject)[i]))  ### specify the covariate matrix
    if(family=="gaussian"){
      var <- switch(corstr,
                    independence=scale*cormax.ind(ncluster), 
                    exchangeable=scale*cormax.exch(ncluster, model.R$geese$alpha),
                    ar1=scale*cormax.ar1(ncluster, model.R$geese$alpha)
      ) }
    if(family=="binomial"){
      var <- switch(corstr,
                    independence=scale*diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster), 
                    exchangeable=scale*diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster)%*%cormax.exch(ncluster, model.R$geese$alpha)%*%diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster),
                    ar1=scale*diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster)%*%cormax.ar1(ncluster, model.R$geese$alpha)%*%diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster)
      )}
    if(family=="poisson"){
      var <- switch(corstr,
                    independence=scale*diag(sqrt(c(exp(covariate%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate%*%beta))),ncluster), 
                    exchangeable=scale*diag(sqrt(c(exp(covariate%*%beta))),ncluster)%*%cormax.exch(ncluster, model.R$geese$alpha)%*%diag(sqrt(c(exp(covariate%*%beta))),ncluster),
                    ar1=scale*diag(sqrt(c(exp(covariate%*%beta))),ncluster)%*%cormax.ar1(ncluster, model.R$geese$alpha)%*%diag(sqrt(c(exp(covariate%*%beta))),ncluster)
      )}
    switch(family,
           gaussian={
           varI <- scale*cormax.ind(ncluster)
           va.ind<-t(covariate)%*%ginv(var)%*%matrix(y_i-mu.R_i)%*%t(matrix(y_i-mu.R_i))%*%ginv(varI)%*%covariate 
           va <- va+va.ind
           },
           binomial={
           varI <-scale*diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate%*%beta)/(1+exp(covariate%*%beta))^2)),ncluster) 
           D<-mat.prod(covariate, exp(covariate%*%beta)/((1+exp(covariate%*%beta))^2))
           va.ind<-t(D)%*%ginv(var)%*%matrix(y_i-mu.R_i)%*%t(matrix(y_i-mu.R_i))%*%ginv(varI)%*%D
           va <- va+va.ind
           },
           poisson={
           varI <-scale*diag(sqrt(c(exp(covariate%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate%*%beta))),ncluster)
           D<-mat.prod(covariate, exp(covariate%*%beta))
           va.ind<-t(D)%*%ginv(var)%*%matrix(y_i-mu.R_i)%*%t(matrix(y_i-mu.R_i))%*%ginv(varI)%*%D
           va <- va+va.ind
           },
           stop("Warnings: Invalid type of outcomes!")
    )
  }	
  trace.R <- sum(diag(va %*% Vnaive))
  # QIC;
  MQIC <- (-2)*quasi.R + 2*trace.R
  MQICu <- (-2)*quasi.R + 2*px
  # output the results;
  return(list(MQIC=MQIC, MQICu=MQICu, Quasi_lik=quasi.R))
}
