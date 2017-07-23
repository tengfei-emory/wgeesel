\name{MQIC.gee}
\alias{MQIC.gee}
\encoding{utf-8}
\title{
MQIC and MQICu for GEE
}
\description{
Calculate MQIC (Modified QIC) and MQIC\eqn{_{u}}{u}. MQIC is an asymptotic unbiased estimator of the risk function based on the independent quasi-likelihood. MQIC and original QIC may have non-negligible effect for model selection, especially when the true correlation structure completely different from the working correlation structure.
}
\usage{
MQIC.gee(object)
}

\arguments{
  \item{object}{a fitted model object of class \code{"wgee"}.}
}

\value{
Return a list of MQIC, MQIC\eqn{_{u}}{u} and Quasi-likelihood.
}
\references{

Pan, W., 2001. Akaike's information criterion in generalized estimating equations. \emph{Biometrics}, 57(1), pp.120-125.

Gosho, M., Hamada, C., and Yoshimura, I. (2011). Modifications of QIC and CIC for Selecting a Working Correlation Structure in the Generalized Estimating Equation Method. \emph{Japanese Journal of Biometrics}, 32(1), 1-12.

Imori, S., 2013. On Properties of QIC in Generalized Estimating Equations. \emph{Hiroshima University}, pp.1-8.
}
\author{
Cong Xu, Zheng Li and Ming Wang
}

\seealso{
\code{\link[geepack]{geeglm}} (\pkg{geepack})
}
\examples{
data(imps)

fit <- wgee(IMPS79 ~ Drug+Sex+Time, data=imps,id=imps$ID,family="gaussian",
             corstr="independence")
MQIC.gee(fit)


data(ohio)

fit2 <- wgee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, 
             family="binomial", corstr="independence")
MQIC.gee(fit2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MQIC}
\keyword{gee}% __ONLY ONE__ keyword per line
