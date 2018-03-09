\name{drgee}
\alias{drgee}
\encoding{utf-8}
\title{
Fit Doubly Robust Weighted Generalized Estimating Equations
}
\description{
 Analyzes longitudinal data with doubly robust augmented GEE approach.
}
\usage{
drgee(model,outcomemodel, data, id, family, corstr, scale = NULL, 
       mismodel = NULL, nameTRT, maxit=200, tol=0.001)
}
\arguments{
  \item{model}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}

\item{outcomemodel}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the ouctome model
}

\item{data}{a data frame containing the variables in the model.}
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
\item{nameTRT}{
name of the variable containing information for the treatment
}
 \item{scale}{
 a numeric variable giving the value to which the scale parameter should be fixed; if NA, the scale parameter is not fixed.
 }
   \item{mismodel}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the missingness model to be fitted. 
}
\item{maxit}{maximum iteration number for Newton-Raphson algorithm.}

   \item{tol}{
the tolorance for the Newton-Raphson algorithm to converge.
}
}

\details{
Analyzes longitudinal data with doubly robust augmented GEE approach. 
}

\value{
An object of type 'CRTgeeDR'

}
\references{

Augmented GEE for improving efficiency and validity of estimation in cluster randomized trials by leveraging cluster-and individual-level covariates -2012-Stephens A., Tchetgen Tchetgen E. and De Gruttola V. : \emph{Stat Med} 31(10)-915-930.

}
\author{
Zheng Li, Cong Xu and Ming Wang
}

\seealso{
https://cran.r-project.org/web/packages/CRTgeeDR/index.html
}
\examples{
data(imps)
fit <- drgee(IMPS79 ~ Drug+Sex+Time, IMPS79 ~ Drug+Sex+Time,
      data=imps, id=imps$ID, family="gaussian",nameTRT="Drug",
       corstr="exchangeable", scale=NULL, mismodel= R ~ Drug+Time)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wgee}
\keyword{MAR} % __ONLY ONE__ keyword per line
