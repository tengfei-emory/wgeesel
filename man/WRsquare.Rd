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
WRsquare.gee(object, weight_mean)
}
\arguments{
  \item{object}{a fitted model object of class \code{"wgee"}.}

\item{weight_mean}{
logical; whether weighted mean of the response should be used for weighted R square.
 }
}

\value{
Return a list of weighted R square.
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

fit <- wgee(Y ~ Drug+Sex+Time, data=imps, id=imps$ID, family="binomial",
            corstr="exchangeable",scale=NULL,mismodel= R ~ Drug+Time)
WRsquare.gee(fit, weight_mean=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wgee}
\keyword{MAR} % __ONLY ONE__ keyword per line
