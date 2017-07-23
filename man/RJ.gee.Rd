\name{RJ.gee}
\alias{RJ.gee}
\encoding{utf-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RJC for GEE
}
\description{
Calculate RJC (Rotnitzky-Jewell information criterion) based on GEE. 
}
\usage{
RJ.gee(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted model object of class \code{"wgee"}.
}
}
\details{
Rotnitzky-Jewell information criterion (RJC) is usually used for working correlation structure selection.
}
\value{
Return the value of the Rotnitzky-Jewell information criterion (RJ). 
}
\references{

Rotnitzky, A. and Jewell, N.P., 1990. Hypothesis testing of regression parameters in semiparametric generalized linear models for cluster correlated data. \emph{Biometrika}, pp.485-497.
}
\author{
Cong Xu, Zheng Li and Ming Wang
}
\seealso{
\code{\link[geepack]{geeglm}} (\pkg{geepack})
}
\examples{
data(ohio)

fit <- wgee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, 
            family="binomial", corstr="exchangeable")
       
RJ.gee(fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RJC}
\keyword{gee}% __ONLY ONE__ keyword per line
