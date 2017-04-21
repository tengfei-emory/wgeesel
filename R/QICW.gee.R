QICW.gee <-
function(model, mismodel, data, id, family, corstr) {

  #########################################################################
  # Arguments: mu.R weight, scale and rho are extracted from macro GEE
  # data     data frame with subject variable from 1:size
  # model    specify the model of interest
  # logit_model  specify the model for dropout
  # corstr   Working correlation structure: "independence", "AR-M", "exchangeable", "unstructured".
  # family    type of the outcomes: gaussian, binomial or poisson
  #### continous case: x*beta;
  #### binomial case: 1/exp(-x*beta)
  #### poisson case: exp(x*beta)
 
  #########################################################################
  # input from wgee()
  # para_est  the estimates of the parameters
  # alpha_drop   the parameter estimates for the dropout model
  # scale  the estimate of variance for gaussian outcomes only based on weighted GEE
  # rho   the estimate of correlation coefficient based on weighted GEE
  # mu.R  the fitted values from the candidate model based on weighted GEE
  # weight   the esimtate of weights based on weighted GEE
  
  #####weight######
  if (length(id) != nrow(data)){
      stop("variable lengths differ (found for '(id)')")}
  m.temp <- model.frame(model, data, na.action='na.pass')
  if (sum(is.na(m.temp)) == 0){
    stop("This data is complete. QICW.gee() is only for data with missingness")}
  
  cluster <-cluster.size(id)$n #number of observations for each subject;
  size <-cluster.size(id)$m #  sample size;
  data$subject <- rep(1:size, cluster) #1:size; create subject variable#
  fit=wgee(model,data,id,family=family,corstr=corstr,mismodel = mismodel)
  data$mu.R=fit$mu_fit
  data$weight=fit$weight
  model=fit$model
  model_drop=fit$mis_fit$formula
  
  para_est=as.vector(fit$beta)
  alpha_drop=as.vector(coef(fit$mis_fit))
  scale=fit$scale
  rho=fit$rho
  init <- model.frame(model, data)
  init$num <- 1:length(init[,1])
  
  ### Get the design matrix;
  m <- model.frame(model, data, na.action='na.pass')
  data$response <- model.response(m, "numeric")
  mat <- as.data.frame(model.matrix(model, m))
  mat$subj <- rep(unique(data$subject), cluster)
     
  beta <- para_est
  px <- length(beta) # number non-redunant columns in design matrix
  
  # Quasi Likelihood;
  index.nmiss <- which(!is.na(data$response))
  switch(family,
         gaussian={model.R <- geeglm(model, id =data$subject, data=data, corstr=corstr, family=gaussian)
         y <- model.R$y
         quasi.R <- -1/2*sum((y-data$mu.R[index.nmiss])^2*data$weight[index.nmiss]) 
         },
         binomial={model.R <- geeglm(model, id =data$subject, data=data, corstr=corstr, family=binomial)        	    
         y <- model.R$y           
         quasi.R <- sum((y*log(data$mu.R[index.nmiss]/(1-data$mu.R[index.nmiss]))+log(1-data$mu.R[index.nmiss]))*data$weight[index.nmiss])
         },
         poisson={model.R <- geeglm(model, id =data$subject, data=data, corstr=corstr, family=poisson)
         y <- model.R$y
         quasi.R <- sum((y*log(data$mu.R[index.nmiss])-data$mu.R[index.nmiss])*data$weight[index.nmiss])
         },
         stop("Warnings: Invalid type of outcomes!")
  )
  
  ##### get the trace term: the penalty for model complexity;
  h.part<-matrix(0, nrow=length(mat[1,])-1, ncol=length(mat[1,])-1) 
  UF<-matrix(0,ncol=length(alpha_drop), nrow=length(mat[1,])-1)
  FF<-matrix(0,ncol=length(alpha_drop), nrow=length(alpha_drop))
  j.list <- k.list <- list()
  for (i in 1:size){
    cluster_alt <- cluster.size(data[!is.na(data$response),]$subject)$n
    ncluster <- cluster_alt[i]
    y<-as.matrix(data[data$subject==i,]$response)
    mu.R<-as.matrix(data[data$subject==i,]$mu.R)
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$subject)[i]))  ### specify the covariate matrix
    switch(family,
           gaussian={
           index <- which(!is.na(data[data$subject==i,]$response))
           var <- switch(corstr,
                         independence=scale*cormax.ind(ncluster), 
                         exchangeable=scale*cormax.exch(ncluster, rho),
                         ar1=scale*cormax.ar1(ncluster, rho)
           )
           if (is.matrix(covariate[index,])==FALSE) {cal.m=t(t(covariate[index,]))
           }else{cal.m=t(covariate[index,])}
           
           h.ind<-cal.m%*%ginv(var)%*%diag(c(data[data$subject==i,]$weight[index]))%*%covariate[index,]
           #h.ind<-t(covariate[index,])%*%ginv(var)%*%covariate[index,]
           h.part<-h.part+h.ind
           ###Calculate the G_i in the formula;
           x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R[index])
           U<-cal.m%*%ginv(var)%*%x
           ### Get the design matrix for dropout model;
           m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
           mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
           exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))      
           
           lam_d <- exp.x/(1+exp.x)
           row.names(m_drop)=seq(1:cluster[i])
           index_d=as.numeric(row.names(na.omit(m_drop)))
           F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ) )
             
           UF<-UF+U%*%t(F)
           FF<-FF+F%*%t(F)                                                    
           },
           binomial={
           index <- which(!is.na(data[data$subject==i,]$response))
           var <- switch(corstr,
                         independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster), 
                         exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster),
                         ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)
           )
           D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta)/((1+exp(covariate[index,]%*%beta))^2))
           h.ind<-t(D)%*%ginv(var)%*%diag(c(data[data$subject==i,]$weight[index]))%*%D
           #h.ind<-t(D)%*%ginv(var)%*%D
           h.part<-h.part+h.ind                
           ###Calculate the G_i in the formula;
           x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R[index])
           U<-t(D)%*%ginv(var)%*%x
           ### Get the design matrix for dropout model;
           m_drop <- model.frame(model_drop, data[data$subject==i,], na.action='na.pass')
           mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
           exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))              
           lam_d <- exp.x/(1+exp.x)
           row.names(m_drop)=seq(1:cluster[i])
           index_d=as.numeric(row.names(na.omit(m_drop)))
            F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ) )
           
           UF<-UF+U%*%t(F)
           FF<-FF+F%*%t(F)    
           },
           poisson={
           index <- which(!is.na(data[data$subject==i,]$response))
           var <- switch(corstr,
                         independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster), 
                         exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster),
                         ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)
           )
           D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta))
           h.ind<-t(D)%*%ginv(var)%*%diag(c(data[data$subject==i,]$weight[index]))%*%D
           h.part<-h.part+h.ind 
           ###Calculate the G_i in the formula;
           x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R[index])
           U<-t(D)%*%ginv(var)%*%x
           ### Get the design matrix for dropout model;
           m_drop <- model.frame(model_drop, data[data$subject==i,], na.action='na.pass')
           mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
           exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))       
           
           
           lam_d <- exp.x/(1+exp.x)
           row.names(m_drop)=seq(1:cluster[i])
           index_d=as.numeric(row.names(na.omit(m_drop)))
          F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ) )
             

           
           UF<-UF+U%*%t(F)
           FF<-FF+F%*%t(F)      
           },
           stop("Warnings: Invalid type of outcomes!")
    )
  }
  G_pre<-UF%*%ginv(FF) ####the first two terms of Gi;
  j.part<-matrix(0,nrow=length(mat[1,])-1,ncol=length(mat[1,])-1)
  k.part<-matrix(0,nrow=length(mat[1,])-1,ncol=length(mat[1,])-1)
  for (i in 1:size){
    cluster_alt <- cluster.size(data[!is.na(data$response),]$subject)$n
    ncluster <- cluster_alt[i]
    y<-as.matrix(data[data$subject==i,]$response)
    mu.R<-as.matrix(data[data$subject==i,]$mu.R)
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$subject)[i]))  ### specify the covariate matrix
    switch(family,
           gaussian={
           index <- which(!is.na(data[data$subject==i,]$response))
           var <- switch(corstr,
                         independence=scale*cormax.ind(ncluster), 
                         exchangeable=scale*cormax.exch(ncluster, rho),
                         ar1=scale*cormax.ar1(ncluster, rho)
           )
           x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R[index])
           ### Get the design matrix for dropout model;
           m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
           mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
           exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
           
           lam_d <- exp.x/(1+exp.x)
           row.names(m_drop)=seq(1:cluster[i])
           index_d=as.numeric(row.names(na.omit(m_drop)))
            F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ) )
           
           G<-G_pre%*%F
           if (is.matrix(covariate[index,])==FALSE) {cal.m=t(t(covariate[index,]))
           }else{cal.m=t(covariate[index,])}
        
           E<-cal.m%*%ginv(var)%*%x
           j.part<-j.part+(E-G)%*%t(E-G)
           varI <- scale*cormax.ind(ncluster)
           #k.part<-k.part+t(covariate[index,])%*%ginv(varI)%*%diag(c(data[data$subject==i,]$weight[index]))%*%covariate[index,]
           k.part<-k.part+cal.m%*%ginv(varI)%*%covariate[index,]
           },
           binomial={
           index <- which(!is.na(data[data$subject==i,]$response))
           var <- switch(corstr,
                         independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster), 
                         exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster),
                         ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)
           )
           D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta)/((1+exp(covariate[index,]%*%beta))^2))              
           x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R[index])
           ### Get the design matrix for dropout model;
           m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
           mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
           exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
           
           lam_d <- exp.x/(1+exp.x)
           row.names(m_drop)=seq(1:cluster[i])
           index_d=as.numeric(row.names(na.omit(m_drop)))
            F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ) )
             
           G<-G_pre%*%F
           E<-t(D)%*%ginv(var)%*%x
           j.part<-j.part+(E-G)%*%t(E-G)
           varI <- scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)
           #k.part<-k.part+t(D)%*%ginv(varI)%*%diag(c(data[data$subject==i,]$weight[index]))%*%D
           k.part<-k.part+t(D)%*%ginv(varI)%*%D
           },
           poisson={
           index <- which(!is.na(data[data$subject==i,]$response))
           var <- switch(corstr,
                         independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster), 
                         exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster),
                         ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)
           )
           D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta))
           x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R[index])
           ### Get the design matrix for dropout model;
           m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
           mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
           exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
           index <- which(!is.na(data[data$subject==i,]$response))
           
           lam_d <- exp.x/(1+exp.x)
           row.names(m_drop)=seq(1:cluster[i])
           index_d=as.numeric(row.names(na.omit(m_drop)))
            F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ) )
            
           G<-G_pre%*%F
           E<-t(D)%*%ginv(var)%*%x
           j.part<-j.part+(E-G)%*%t(E-G)
           varI <- scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)
           #k.part<-k.part+t(D)%*%ginv(varI)%*%diag(c(data[data$subject==i,]$weight[index]))%*%D
           k.part<-k.part+t(D)%*%ginv(varI)%*%D
           },
           stop("Warnings: Invalid type of outcomes!")
    )
  }
  trace.R <- sum(diag(ginv(h.part)%*%j.part%*%ginv(h.part)%*%k.part))
  #trace.R <- sum(diag(k.part%*%ginv(h.part)%*%j.part%*%ginv(h.part)))
  #WQIC;
  WQIC <- (-2)*quasi.R + 2*trace.R
  WQICu <- (-2)*quasi.R + 2*px
  # output the results;
  return(list(QICWr=WQIC, QICWp=WQICu, Quasi_lik=quasi.R))
}
