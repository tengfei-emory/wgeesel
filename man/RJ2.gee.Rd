\name{RJ2.gee}
\alias{RJ2.gee}
\encoding{utf-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Corrected RJC for GEE
}
\description{
Calculate corrected RJC (Rotnitzky-Jewell information criterion) based on GEE with a modified robust variance estimator. 
}
\usage{
RJ2.gee(model, data, id, family, corstr)
}
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
\value{
Return the value of the corrected Rotnitzky-Jewell information criterion (RJC). 
}
\references{

Rotnitzky, A. and Jewell, N.P., 1990. Hypothesis testing of regression parameters in semiparametric generalized linear models for cluster correlated data. \emph{Biometrika}, pp.485-497.

Wang, M. and Long, Q., 2011. Modified robust variance estimator for generalized estimating equations with improved small-sample performance. \emph{Statistics in Medicine}, 30(11), pp. 1278-1291.

}
\author{
Cong Xu, Zheng Li and Ming Wang
}

\note{
\code{\link{RJ2.gee}} can only handle balanced data (data with dropout missingness). Two assumptions should be satisfied. (A1) The conditional variance of \eqn{Y_{ij}} given \eqn{X_{ij}} is correctly specified; (A2) A common correlation structure, \eqn{Rc}, exists across all subjects. If there is missingness, one can group the subjects by the cluter size of the response variable. And, calculate the modified robust variance in each group to get the pooled estimate of the variance. 

}
\seealso{
\code{\link[geepack]{geeglm}} (\pkg{geepack}), \code{\link{RJ.gee}}
}
\examples{
data(ohio)

RJ2.gee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, family="binomial", 
       corstr="exchangeable")
       
RJ2.gee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, family="binomial", 
       corstr="ar1")
}
\keyword{RJC}
\keyword{gee}
