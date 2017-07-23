\name{QIC.gee}
\alias{QIC.gee}
\encoding{utf-8}

\title{
QIC and QICu for GEE
}
\description{
Calculate quasi-likelihood under the independence model criterion (QIC) and QIC\eqn{_{u}}{u} based on GEE.
}
\usage{
QIC.gee(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model object of class \code{"wgee"}.
}
}

\value{
Return a vector of QIC, QIC\eqn{_{u}}{u} and Quasi-likelihood.
}
\references{

Liang, K.Y. and Zeger, S.L., 1986. Longitudinal data analysis using generalized linear models. \emph{Biometrika}, pp.13-22.

Pan, W., 2001. Akaike's information criterion in generalized estimating equations. \emph{Biometrics}, 57(1), pp.120-125.

Prentice, R.L. and Zhao, L.P., 1991. Estimating equations for parameters in means and covariances of multivariate discrete and continuous responses. \emph{Biometrics}, pp.825-839.
}
\author{
Cong Xu, Zheng Li and Ming Wang
}
\note{
QIC can be used to select the best correlation structure and the best fitting model in GEE analyses. The GEE is fitted by \code{\link[geepack]{geeglm}} (\pkg{geepack}). QIC\eqn{_{u}}{u} is a simplified version of QIC, which can not be applied to select the optimal working correlation structure. \code{\link[geepack]{geeglm}} (\pkg{geepack}) only works for complete data. Thus if there are NA’s in data, the missing values are automatically removed by \code{\link[stats]{na.omit}}.
}

\seealso{
\code{\link[geepack]{geeglm}} (\pkg{geepack}). \pkg{MuMIn} also provides QIC value.
}

\examples{
data(imps)
fit <- wgee(Y ~ Drug+Sex+Time, data=imps, id=imps$ID, family="binomial",
            corstr="exchangeable")

QIC.gee(fit)

data(seizure)
###reshapre the seizure data to "long" format
seiz.long <- reshape(seizure, varying=list(c("base","y1", "y2", "y3", "y4")),
                      v.names="y", times=0:4, direction="long")
seiz.long <- seiz.long[order(seiz.long$id, seiz.long$time),]

fit <- wgee(y ~ age + trt + time, data=seiz.long, id=seiz.long$id, 
            family="poisson", corstr="independence")

QIC.gee(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QIC}
\keyword{gee}% __ONLY ONE__ keyword per line
