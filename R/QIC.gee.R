QIC.gee <-
function(object) {

  model=object$model
  data=object$data
  id=object$id
  corstr=object$corstr
  family=object$family
  
  ### data information;
  if (length(id) != nrow(data)){
    stop("variable lengths differ (found for '(id)')")}
  m <- model.frame(model, data, na.action='na.pass')
  cluster<-cluster.size(id)$n #number of observations for each subject;
  size<-cluster.size(id)$m # sample size;
  
  subject <- rep(1:size, cluster) #1:size;
  ### Get the design matrix;
  data.temp <- cbind(subject,m)
  ### get rid of the missing data;
  data <- na.omit(data.temp)
  # Quasi Likelihood;
  switch(family,
         gaussian={model.ind <- geeglm(model, id =subject, data=data, corstr="independence", family=gaussian)
         ### alternative correlation structure;
         model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=gaussian)
         #model.R <- gee(model, id =subject, data=data, corstr=corstr, family=gaussian)
         mu.R <- model.R$fitted.values
         y <- model.R$y
         scale <- as.numeric(summary(model.R)$geese$scale[1])
         quasi.R <- -1/2*sum((y-mu.R)^2/scale)
         },
         binomial={model.ind <- geeglm(model, id =subject, data=data, corstr="independence", family=binomial)
         ### alternative correlation structure;
         model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=binomial)
         mu.R <- model.R$fitted.values
         y <- model.R$y
         quasi.R <- sum(y*log(mu.R/(1-mu.R))+log(1-mu.R))
         },
         poisson={model.ind <- geeglm(model, id =subject, data=data, corstr="independence", family=poisson)
         ### alternative correlation structure;
         model.R <- geeglm(model, id =subject, data=data, corstr=corstr, family=poisson)
         mu.R <- model.R$fitted.values
         y <- model.R$y
         quasi.R <- sum(y*log(mu.R)-mu.R)
         },
         stop("Warnings: Invalid type of outcomes!")
  )
  # Trace Term (penalty for model complexity);
  AIinverse <- solve(model.ind$geese$vbeta.naiv) 
  # Omega-hat(I) via Moore-Penrose generalized inverse of a matrix;
  # Alt: AIinverse <- solve(model.ind$vbeta.naiv) # solve via identity;
  Vr <- model.R$geese$vbeta
  trace.R <- sum(diag(AIinverse %*% Vr))
  px <- dim(Vr)[1] # number non-redunant columns in design matrix
  # QIC;
  QIC <- (-2)*quasi.R + 2*trace.R
  QICu <- (-2)*quasi.R + 2*px
  # output the results;
  out <- data.frame(QIC=QIC, QICu=QICu, Quasi_lik=quasi.R)
  return(round(out,1))
  #return(list(QIC=QIC, QICu=QICu, Quasi_lik=quasi.R))
}
