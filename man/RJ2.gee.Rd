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
RJ2.gee(object)
}
\arguments{
  \item{object}{
a fitted model object of class \code{"wgee"}.
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

fit <- wgee(resp ~ age + smoke+age:smoke, data=ohio, id=ohio$id, 
            family="binomial",corstr="exchangeable")
RJ2.gee(fit)       

}
\keyword{RJC}
\keyword{gee}
