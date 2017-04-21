MLIC.gee <-
function(model, model_full, mismodel, data, id, family,corstr) {
  #########################################################################
  # Arguments: mu.R, mu.R0, weight, scale and rho are extracted from macro GEE
  # model    specify the model of interest
  # model_full  specify the full model
  # model_drop  specify the model for dropout
  # data     data frame with subject variable from 1:size
  # id       specify the id
  # family   type of the outcomes: gaussian, binomial or poisson
  #### continous case: x*beta;
  #### binomial case: 1/exp(-x*beta)
  #### poisson case: exp(x*beta)
  # corstr   Working correlation structure: "independence", "AR-M", "exchangeable", "unstructured".
  # complete   the indictor to show if the data is complete or not. 
  #########################################################################
  # para_est  the estimates of the parameters
  # alpha_drop   the parameter estimates for the dropout model
  # scale  the estimate of variance for gaussian outcomes only based on weighted GEE
  # mu.R  the fitted values from the candidate model based on weighted GEE
  # mu.R0   the fitted values from the full model based on weighted GEE
  # weight   the esimtate of weights based on weighted GEE (this should be 1/wij)
  
  #input mu.R, mu.R0, weight, scale, rho, para_est, alpha_drop from wgee###
  m.temp <- model.frame(model_full, data, na.action='na.pass')
  if (sum(is.na(m.temp)) == 0){
    stop("This data is complete. MLIC.gee() is only for data with missingness")}
  complete=F
  if (length(id) != nrow(data)){
      stop("variable lengths differ (found for '(id)')")}
 
  fit=wgee(model,data,id,family=family,corstr=corstr,mismodel=mismodel)
  weight=fit$weight
  full_fit=wgee(model_full,data,id,family=family,corstr=corstr,mismodel=mismodel)
  data$mu.R=fit$mu_fit
  data$mu.R0=full_fit$mu_fit
  data$weight=fit$weight
  model=fit$model
  model_full=full_fit$model
  model_drop=fit$mis_fit$formula
  
  para_est=as.vector(fit$beta)
  alpha_drop=as.vector(coef(fit$mis_fit))
  scale=fit$scale
  rho=fit$rho
  init <- model.frame(model, data)
  init$num <- 1:length(init[,1])
  
  cluster <-cluster.size(id)$n #number of observations for each subject;
  size <-cluster.size(id)$m #  sample size;
  data$subject <- rep(1:size, cluster) #1:size; create subject variable#
  
  ### Get the design matrix;
  m <- model.frame(model, data, na.action='na.pass')
  data$response <- model.response(m, "numeric")
  mat <- as.data.frame(model.matrix(model, m))
  mat$subj <- rep(unique(data$subject), cluster)
  complete=F
  if (complete==T) {
    ### Specify the type of the outcomes; Here requires "geepack" package;
    switch(family,
           gaussian={model.R <- geeglm(model, id=data$subject, data=data, corstr=corstr, family=gaussian)
           model.full <- geeglm(model_full, id=data$subject, data=data, corstr=corstr, family=gaussian)
           data$mu.R <- model.R$fitted.values
           data$mu.R0 <- model.full$fitted.values
           y <- model.R$y
           beta <- as.numeric(model.R$coefficients)
           },
           binomial={model.R <- geeglm(model, id=data$subject, data=data, corstr=corstr, family=binomial)
           model.full <- geeglm(model_full, id =data$subject, data=data, corstr=corstr, family=binomial)
           data$mu.R <- model.R$fitted.values
           data$mu.R0 <- model.full$fitted.values
           y <- model.R$y              
           beta <- as.numeric(model.R$coefficients)
           },
           poisson={model.R <- geeglm(model, id =data$subject, data=data, corstr=corstr, family=poisson)
           model.full <- geeglm(model_full, id =data$subject, data=data, corstr=corstr, family=poisson)
           data$mu.R <- model.R$fitted.values
           data$mu.R0 <- model.full$fitted.values
           y <- model.R$y              
           beta <- as.numeric(model.R$coefficients)
           },
           stop("Warnings: Invalid type of outcomes!")
    )
    # The quasi likelihood part;
    quasi.R <- (y-data$mu.R)%*%matrix(y-data$mu.R) 
    
    # Trace Term (penalty for model complexity);
    h.part <- matrix(0, nrow=length(beta), ncol=length(beta)) #p by p matrix
    j.part <- matrix(0,nrow=length(beta),ncol=length(beta)) #p by p matrix
    
    for (i in 1:size){
      ncluster <- cluster[i]
      y<-as.matrix(data[data$subject==i,]$response)
      mu.R<-as.matrix(data[data$subject==i,]$mu.R)
      mu.R0<-as.matrix(data[data$subject==i,]$mu.R0)
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
             gaussian={h.ind<-t(covariate)%*%ginv(var)%*%covariate
             h.part<-h.part+h.ind
             j.ind<-t(covariate)%*%ginv(var)%*%matrix(y-mu.R0)%*%t(matrix(y-mu.R0))%*%covariate 
             j.list[[i]] <- j.ind   
             j.part <- j.part+j.ind
             },
             binomial={D<-mat.prod(covariate, exp(covariate%*%beta)/((1+exp(covariate%*%beta))^2))
             h.ind <-t(D)%*%ginv(var)%*%D
             h.part<-h.part+h.ind
             j.ind <-t(D)%*%ginv(var)%*%matrix(y-mu.R0) %*%t(matrix(y-mu.R0))%*%D 
             j.part <- j.part+j.ind
             },
             poisson={D<-mat.prod(covariate, exp(covariate%*%beta))
             h.ind<-t(D)%*%ginv(var)%*%D
             h.part<-h.part+h.ind
             j.ind<-t(D)%*%ginv(var)%*%matrix(y-mu.R0) %*%t(matrix(y-mu.R0))%*%D 
             j.part<-j.part+j.ind
             },
             stop("Warnings: Invalid type of outcomes!")
      )
    }	
    # missing information criteria
    mlic <- formatC(quasi.R+2*sum(diag(ginv(h.part)%*%j.part)), digits=3, format="f") 
  } 
  
  if (complete==F) {
    beta <- para_est
    h.part<-matrix(0, nrow=length(mat[1,])-1, ncol=length(mat[1,])-1) 
    UF<-matrix(0,ncol=length(alpha_drop), nrow=length(mat[1,])-1)
    FF<-matrix(0,ncol=length(alpha_drop), nrow=length(alpha_drop))
    j.list <- k.list <- list()
    for (i in 1:size){
      cluster_alt <- cluster.size(data[!is.na(data$response),]$subject)$n
      ncluster <- cluster_alt[i]
      y<-as.matrix(data[data$subject==i,]$response)
      mu.R<-as.matrix(data[data$subject==i,]$mu.R)
      mu.R0<-as.matrix(data[data$subject==i,]$mu.R0)
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$subject)[i]))  ### specify the covariate matrix
      
      switch(family,
             gaussian={
             index <- which(!is.na(data[data$subject==i,]$response))
             var <- switch(corstr,
                           independence=scale*cormax.ind(ncluster), 
                           exchangeable=scale*cormax.exch(ncluster, rho),
                           ar1=scale*cormax.ar1(ncluster, rho)
                         )
             ##covariate[index,]is not a matrix. t(covariate) is row vector. we need col vector.
             if (is.matrix(covariate[index,])==FALSE) {cal.m=t(t(covariate[index,]))
             }else{cal.m=t(covariate[index,])}
             h.ind<-cal.m%*%ginv(var)%*%diag(c(data[data$subject==i,]$weight[index]))%*%covariate[index,]
             h.part<-h.part+h.ind
             ###Calculate the G_i in the formula;
             x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R0[index])
             U<-cal.m%*%ginv(var)%*%x
             ### Get the design matrix for dropout model;
             m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
             mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
             exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
             
             lam_d <- exp.x/(1+exp.x)
             row.names(m_drop)=seq(1:cluster[i])
             index_d=as.numeric(row.names(m_drop))[-1]
             F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ),na.rm=T )
             
             UF<-UF+U%*%t(F)
             FF<-FF+F%*%t(F)                                                    
             }  ,
        
             binomial={
             index <- which(!is.na(data[data$subject==i,]$response))
             var <- switch(corstr,
                           independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster), 
                           exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster),
                           ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)
             )
             D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta)/((1+exp(covariate[index,]%*%beta))^2))
             h.ind<-t(D)%*%ginv(var)%*%diag(c(data[data$subject==i,]$weight[index]))%*%D
             h.part<-h.part+h.ind                
             ###Calculate the G_i in the formula;
             x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R0[index])
             U<-t(D)%*%ginv(var)%*%x
             ### Get the design matrix for dropout model;

             m_drop <- model.frame(model_drop, data[data$subject==i,], na.action='na.pass')
             mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
             exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))              
             lam_d <- exp.x/(1+exp.x)
             row.names(m_drop)=seq(1:cluster[i])
             index_d=as.numeric(row.names(m_drop))[-1]
             F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,],na.rm=T) )
             UF<-UF+U%*%t(F)
             FF<-FF+F%*%t(F)  
             },
             poisson={
             index <- which(!is.na(data[data$subject==i,]$response))
             var <- switch(corstr,
                           independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster), 
                           exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(c(sqrt(exp(covariate[index,]%*%beta))),ncluster),
                           ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)
             )
             D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta))
             h.ind<-t(D)%*%ginv(var)%*%diag(c(data[data$subject==i,]$weight[index]))%*%D
             h.part<-h.part+h.ind 
             ###Calculate the G_i in the formula;
             x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R0[index])
             U<-t(D)%*%ginv(var)%*%x
             ### Get the design matrix for dropout model;
             m_drop <- model.frame(model_drop, data[data$subject==i,], na.action='na.pass')
             mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
             exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))      
             
             lam_d <- exp.x/(1+exp.x)
             row.names(m_drop)=seq(1:cluster[i])
             index_d=as.numeric(row.names(m_drop))[-1]
         
             F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,] ,na.rm=T) )
             
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
      mu.R0<-as.matrix(data[data$subject==i,]$mu.R0)
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$subject)[i]))  ### specify the covariate matrix
      ### get the Omega matrix for mlicc;
      if (ncluster==1) {Omega= 1
      }else{Omega <- matrix(1,nrow=ncluster, ncol=ncluster)
      for (index in 2:ncluster){
        Omega[index:ncluster, index:ncluster] <- matrix(1/data[data$subject==i,]$weight[index], nrow=ncluster-index+1, ncol=ncluster-index+1)
      }
      }
      switch(family,
             gaussian={
             index <- which(!is.na(data[data$subject==i,]$response))
             var <- switch(corstr,
                           independence=scale*cormax.ind(ncluster), 
                           exchangeable=scale*cormax.exch(ncluster, rho),
                           ar1=scale*cormax.ar1(ncluster, rho)
             )
             x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R0[index])
             ### Get the design matrix for dropout model;
             m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
             mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
             exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
             
             lam_d <- exp.x/(1+exp.x)
             row.names(m_drop)=seq(1:cluster[i])
             index_d=as.numeric(row.names(m_drop))[-1]
             F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,],na.rm=T ) )
           
             G<-G_pre%*%F
             if (is.matrix(covariate[index,])==FALSE) {cal.m=t(t(covariate[index,]))
             }else{cal.m=t(covariate[index,])}
             first<-cal.m%*%ginv(var)%*%(x%*%t(x))
             second<-G%*%t(x)
             j.part<-j.part+(first-second)%*%covariate[index,]
             k.part<-k.part+cal.m%*%ginv(var)%*%(Omega*(x%*%t(x)))%*%covariate[index,]
             },
             
             binomial={
             index <- which(!is.na(data[data$subject==i,]$response))
             var <- switch(corstr,
                           independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster), 
                           exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster),
                           ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta)/(1+exp(covariate[index,]%*%beta))^2)),ncluster)
             )
             D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta)/((1+exp(covariate[index,]%*%beta))^2))              
             x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R0[index])
             ### Get the design matrix for dropout model;
             m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
             mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
             exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
          
             lam_d <- exp.x/(1+exp.x)
             row.names(m_drop)=seq(1:cluster[i])
             index_d=as.numeric(row.names(m_drop))[-1]
             F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,],na.rm=T ) )
             
             G<-G_pre%*%F
             first<-t(D)%*%ginv(var)%*%(x%*%t(x))
             second<-G%*%t(x)
             j.part<-j.part+(first-second)%*%D 
             j.list[[i]] <- j.part
             k.part<-k.part+t(D)%*%ginv(var)%*%(Omega*(x%*%t(x)))%*%D    
             k.list[[i]] <- k.part
             },
             poisson={
             index <- which(!is.na(data[data$subject==i,]$response))
             var <- switch(corstr,
                           independence=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ind(ncluster)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster), 
                           exchangeable=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.exch(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster),
                           ar1=scale*diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)%*%cormax.ar1(ncluster, rho)%*%diag(sqrt(c(exp(covariate[index,]%*%beta))),ncluster)
             )
             D<-mat.prod(covariate[index,], exp(covariate[index,]%*%beta))
             x<-diag(c(data[data$subject==i,]$weight[index]))%*%matrix(y[index]-mu.R0[index])
             ### Get the design matrix for dropout model;
             m_drop <- model.frame(model_drop, data[data$subject==i,],na.action='na.pass')
             mat_drop <- as.data.frame(model.matrix(model_drop, m_drop))
             exp.x <- exp(as.matrix(mat_drop)%*%as.matrix(alpha_drop))
      
             lam_d <- exp.x/(1+exp.x)
             row.names(m_drop)=seq(1:cluster[i])
             index_d=as.numeric(row.names(m_drop))[-1]
             F= as.matrix(colSums( (m_drop[index_d-1,][,1] - m_drop[-index_d[length(index_d)],][,1]*c(lam_d[index_d-1]) )*
                                     mat_drop[index_d-1,],na.rm=T ) )
             
             G<-G_pre%*%F
             first<-t(D)%*%ginv(var)%*%(x%*%t(x))
             second<-G%*%t(x)
             j.part<-j.part+(first-second)%*%D 
             k.part<-k.part+t(D)%*%ginv(var)%*%(Omega*(x%*%t(x)))%*%D        
             },
             stop("Warnings: Invalid type of outcomes!")
      )
    }
    # The quasi likelihood part;
    index.nmiss <- which(!is.na(data$response))
    quasi.R <- (data$response[index.nmiss]-data$mu.R[index.nmiss])%*%diag(c(data$weight[index.nmiss]))%*%matrix(data$response[index.nmiss]-data$mu.R[index.nmiss]) 
    # missing information criteria: mlic
    mlic <- formatC(quasi.R+2*sum(diag(ginv(h.part)%*%j.part)), digits=3, format="f")
    # missing informaiton criteria for correlation structure: mlicc
    mlicc <- formatC(quasi.R+2*sum(diag(ginv(h.part)%*%k.part)), digits=3, format="f")
  }
  # output the results;
  return(c(MLIC=mlic, MLICc=mlicc, Quasi_lik=quasi.R))
}
