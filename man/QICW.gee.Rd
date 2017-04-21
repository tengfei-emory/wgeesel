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
QICW.gee(model, mismodel, data, id, family, corstr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
}
  \item{mismodel}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the missingness model to be fitted. 
}
  \item{data}{
a data frame containing the variables in the model.
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
Return a list of QICW\eqn{_{r}}{r}, QICW\eqn{_{p}}{p} and Quasi-likelihood.
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

### variable selection by QICWr
QICW.gee(IMPS79 ~ Drug+Sex+Time, mismodel= R ~ Drug+Time, data=imps, id=imps$ID,
         family="gaussian", corstr="exchangeable")

QICW.gee(IMPS79 ~ Drug+Sex+Time+Time:Sex+Time:Drug, mismodel= R ~ Drug+Time,
         data=imps, id=imps$ID, family="gaussian", corstr="exchangeable")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QICW}
\keyword{wgee}% __ONLY ONE__ keyword per line
