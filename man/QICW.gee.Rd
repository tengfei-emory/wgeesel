\name{QICW.gee}
\alias{QICW.gee}
\encoding{utf-8}

\title{
QICWr and QICWp for WGEE
}
\description{
Calculate the QICW\eqn{_{r}}{r} and QICW\eqn{_{p}}{p} (an information criterion based on the weighted quasi-likelihood function) for selection of mean model and correlation structure based on the WGEE.
}
\usage{
QICW.gee(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model object of class \code{"wgee"}.}
}

\value{
Return a data frame of QICW\eqn{_{r}}{r}, QICW\eqn{_{p}}{p} and Wquasi_lik.
}
\references{
Gosho, M., 2015. Model selection in the weighted generalized estimating equations for longitudinal data with dropout. \emph{Biometrical Journal}.

Platt, R.W., Brookhart, M.A., Cole, S.R., Westreich, D. and Schisterman, E.F., 2013. An information criterion for marginal structural models. \emph{Statistics in medicine}, 32(8), pp.1383-1393.

Robins, J.M., Rotnitzky, A. and Zhao, L.P., 1995. Analysis of semiparametric regression models for repeated outcomes in the presence of missing data. \emph{Journal of the American Statistical Association}, 90(429), pp.106-121.
}

\author{
Cong Xu, Zheng Li and Ming Wang
}

\note{

QICW\eqn{_{r}}{r} can be used for variable selection and for selecting the correlation structure in WGEE analyses. QICW\eqn{_{p}}{p} is a simplified version of QICW\eqn{_{r}}{r}, which can not be applied to select the optimal working correlation structure in WGEE.

} 

\seealso{
\code{\link{wgee}}
}

\examples{
data(imps)

### variable selection by QICWr, not rum###
#fit <- wgee(IMPS79 ~ Drug+Sex+Time, mismodel= R ~ Drug+Time, data=imps,
##           id=imps$ID, family="gaussian", corstr="exchangeable")
##QICW.gee(fit)

#fit <- wgee(IMPS79 ~ Drug+Sex+Time+Time:Sex+Time:Drug, mismodel= R ~ Drug+Time,
#        data=imps, id=imps$ID, family="gaussian", corstr="exchangeable")
##QICW.gee(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QICW}
\keyword{wgee}% __ONLY ONE__ keyword per line
