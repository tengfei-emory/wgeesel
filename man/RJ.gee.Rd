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
RJ.gee(model, data, id, family, corstr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
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

RJ.gee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, family="binomial", 
       corstr="exchangeable")
       
RJ.gee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, family="binomial", 
       corstr="ar1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RJC}
\keyword{gee}% __ONLY ONE__ keyword per line
