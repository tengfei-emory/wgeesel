\name{MLIC.gee}
\alias{MLIC.gee}
\encoding{utf-8}
\title{
MLIC and MLICC for Weighted GEE
}
\description{
Calculate the MLIC (missing longitudinal information
criterion) for selection of mean model, and the MLICC (missing longitudinal information correlation
criterion) for selection of working correlation structure, based on the expected quadratic loss and the WGEE.
}
\usage{
MLIC.gee(model, model_full, mismodel, data, id, family, corstr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
}
  \item{model_full}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the largest candidate model under consideration to be fitted. 
}
  \item{mismodel}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the missingness model to be fitted. 
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
"independence", "exchangeable", "ar1".
}
}
 
\value{
Return a data frame of MLIC, MLICC and Wquad_loss.
}

\author{
Cong Xu, Zheng Li and Ming Wang
}

\references{
Robins, J.M., Rotnitzky, A. and Zhao, L.P., 1995. Analysis of semiparametric regression models for repeated outcomes in the presence of missing data. \emph{Journal of the American Statistical Association}, 90(429), pp.106-121.

Shen, C.W. and Chen, Y.H., 2012. Model selection for generalized estimating equations accommodating dropout missingness. \emph{Biometrics}, 68(4), pp.1046-1054.

Shen, C.W. and Chen, Y.H., 2013. Model selection of generalized estimating equations with multiply imputed longitudinal data. \emph{Biometrical Journal}, 55(6), pp.899-911.
}
\note{
MLIC and MLICC model selection criterion for longitudinal data criterion with dropouts or monotone missingness under the assumption of MAR.
}

\seealso{
\code{\link{wgee}}
}

\examples{
data(imps)

###create the lag variable of Y
ylag.list <- lapply(unique(imps$ID),function(x){idx=which(imps$ID==x)
ylag=c(NA,imps$Y[idx][1:(length(idx)-1)])})
imps$ylag <- unlist(ylag.list)

###fit WGEE
fit <- wgee(Y ~ Drug+Sex+Time, data=imps, id=imps$ID, family="binomial",
            corstr="exchangeable", scale=NULL, mismodel= R ~ Drug+Time)

### variable selection by MLIC
#MLIC.gee(Y ~ Drug+Time, Y ~ Drug+Sex+Time+Time:Drug+Time:Sex,
#         mismodel= R ~ Drug+Time+ylag, data=imps, id=imps$ID,family="binomial",
#         corstr="ar1")

#MLIC.gee(Y ~ Drug+Sex+Time, Y ~ Drug+Sex+Drug:Sex+Time:Drug+Time:Sex, 
#         mismodel= R ~ Drug+Time+ylag, data=imps, id=imps$ID, family="binomial",
#         corstr="ar1")

#MLIC.gee(Y ~ Drug+Sex+Time+Drug:Sex, Y ~ Drug+Sex+Drug:Sex+Time:Drug+Time:Sex,
#         mismodel= R ~ Drug+Time+ylag, data=imps, id=imps$ID, family="binomial",
#         corstr="ar1")

}

\keyword{MLIC}
\keyword{wgee}% __ONLY ONE__ keyword per line
