\name{wgee}
\alias{wgee}
\encoding{utf-8}
\title{
Fit Weighted Generalized Estimating Equations (WGEE)
}
\description{
\code{\link{wgee}} fits weighted generalized estimating equations (WGEE) with Newton Raphson algorithm. \code{\link{wgee}} has a syntax similar to \code{\link[stats]{glm}} and returns an object similar to a \code{\link[stats]{glm}} object. 
}
\usage{
wgee(model, data, id, family, corstr, scale = NULL, mismodel = NULL, maxit=200, tol=0.001)
}
\arguments{
  \item{model}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{data}{a data frame containing the variables in the model.
}
 \item{id}{
a vector which identifies the clusters. The length of "id" should be the same as
the number of observations. Data are assumed to be sorted so that observations
on a cluster are contiguous rows for all entities in the formula. 
}
  \item{family}{
a description of the error distribution and link function to be used in the model. This is a character string naming a family function. The following are permitted: "gaussian", "binomial", "poisson".
}
  \item{corstr}{
a character string specifies the working correlation structure. The following are permitted:
"independence", "exchangeable", "ar1","unstructured".
}
 \item{scale}{
 a numeric variable giving the value to which the scale parameter should be fixed; if NA, the scale parameter is not fixed.
 }
   \item{mismodel}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the missingness model to be fitted. 
}

   \item{maxit}{
maximum iteration number for Newton-Raphson algorithm.
}

   \item{tol}{
the tolorance for the Newton-Raphson algorithm to converge.
}
}

\details{
\code{\link{wgee}} analyzes longitudinal data with missing values by weighted genralized estimating equations (WGEE), proposed by Robins, Totnizky and Zhao (1995). WGEE can handle missing at random problem. The standard error of the estimates are calculated as described in (Fitzmaurice, Laird, and Ware, 2011) and Preisser, Lohman, and Rathouz (2002).
}

\value{
\item{beta}{ covariate effect estimates}
\item{var}{ variance covariances estimates for beta}
\item{w_r_square}{ weighted R square for continuous data}
\item{mu_fit}{ fitted values of response}
\item{scale}{scale estimates}
\item{rho}{rho estimates}
\item{weight}{The weight of response y}
\item{model}{WGEE model structure}
\item{x}{covariates in WGEE}
\item{y}{response in WGEE}
\item{mis_fit}{esimates of the missingness model}
\item{call}{the function to be called}
\item{id}{as input}
\item{data}{as input}
\item{family}{as input}
\item{corstr}{as input}

}
\references{

Fitzmaurice, G.M., Laird, N.M. and Ware, J.H., 2012. \emph{Applied longitudinal analysis (Vol. 998)}. John Wiley & Sons.

Liang, K.Y. and Zeger, S.L., 1986. Longitudinal data analysis using generalized linear models. \emph{Biometrika}, pp.13-22.

Preisser, J.S., Lohman, K.K. and Rathouz, P.J., 2002. Performance of weighted estimating equations for longitudinal binary data with drop-outs missing at random. \emph{Statistics in medicine}, 21(20), pp.3035-3054.

Robins, J.M., Rotnitzky, A. and Zhao, L.P., 1995. Analysis of semiparametric regression models for repeated outcomes in the presence of missing data. \emph{Journal of the American Statistical Association}, 90(429), pp.106-121.

Rubin, D.B., 1976. Inference and missing data. \emph{Biometrika}, pp.581-592.

}
\author{
Zheng Li, Cong Xu and Ming Wang
}
\seealso{
\code{\link[geepack]{geeglm}} (\pkg{geepack})
}

\examples{
####Example1
data(imps)

fit <- wgee(IMPS79 ~ Drug+Sex+Time, data=imps, id=imps$ID, family="gaussian",
            corstr="exchangeable", scale=NULL, mismodel= R ~ Drug+Time)

#####Example2
data(seizure)

###reshapre the seizure data to "long" format
seiz.long <- reshape(seizure,
                     varying=list(c("base","y1", "y2", "y3", "y4")),
                     v.names="y", times=0:4, direction="long")
seiz.long <- seiz.long[order(seiz.long$id, seiz.long$time),]

###create missing value for seiz.long dataset
set.seed(12345)
obs <- exp(9+seiz.long$age*(-0.2))/(1+exp(9+seiz.long$age*(-0.2)))
R <- lapply(unique(seiz.long$id),function(x){
  idx=which(seiz.long$id==x)
  r=c()
  r[1]=1
  for(j in 2:length(idx)){
    if(r[j-1]==1){
      r[j]=rbinom(1,1,obs[idx[j]])
    }
    else r[j]=0
  }
  return(r)
})
remove_id <- which(sapply(R,sum)==1)
remove_idx <- which(seiz.long$id \%in\% remove_id==1)
seiz.long <- cbind(seiz.long,R=unlist(R))[-remove_idx,]
seiz.long$y_mis <- ifelse(seiz.long$R,seiz.long$y,NA)

###fit WGEE, not run ###
fit <- wgee(y_mis ~ age + trt + time, data=seiz.long, id=seiz.long$id,family="poisson",
            corstr="exchangeable",scale=NULL, mismodel= R ~ age)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wgee}
\keyword{MAR} % __ONLY ONE__ keyword per line
