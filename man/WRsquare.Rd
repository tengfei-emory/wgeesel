\name{WRsquare.gee}
\alias{WRsquare.gee}
\encoding{utf-8}
\title{
Weighted R Square for WGEE
}
\description{
Calculate the weighted \eqn{R^2} (missing longitudinal information criterion) for selection of mean model.
}
\usage{
WRsquare.gee(model, data, id, family, corstr, scale, mismodel, weight_mean)
}
\arguments{
  \item{model}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
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
"independence", "exchangeable", "ar1".
}
 \item{scale}{
 a numeric variable giving the value to which the scale parameter should be fixed; if NA, the scale parameter is not fixed.
 }
\item{mismodel}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the missingness model to be fitted. 
}

\item{weight_mean}{
logical; whether weighted mean of the response should be used for weighted R square.
 }
}

\value{
Return a list of weighted R square
}
\references{
Nelder, J.A. and Baker, R.J., 1972. \emph{Generalized linear models}. Encyclopedia of statistical sciences.
}
\author{
Zheng Li, Cong Xu and Ming Wang
}
\seealso{
\code{\link{wgee}}
}

\examples{
data(imps)

WRsquare.gee(Y ~ Drug+Sex+Time, data=imps, id=imps$ID, family="binomial",
             corstr="exchangeable",scale=NULL,mismodel= R ~ Drug+Time, 
             weight_mean=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wgee}
\keyword{MAR} % __ONLY ONE__ keyword per line
