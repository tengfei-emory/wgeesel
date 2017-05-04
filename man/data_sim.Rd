\name{data_sim}
\alias{data_sim}
\encoding{utf-8}
\title{
Simulate longitudinal data
}
\description{
The function generate correlated normal, Bernoulli or Poisson longitudinal data.
}
\usage{
data_sim(id, rho, phi, x, beta, x_mis, para, corstr, family, lag_level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
subject id
}
  \item{rho}{
with-in cluster correlation.
}
  \item{phi}{
scale parameter in the variance covariance matrix.
}
  \item{x}{
covariate associated with the response.
}
  \item{beta}{
coefficients associated with x.
}
  \item{x_mis}{
covariates associated with missing model.
}
  \item{para}{
coefficients associated with x_mis.
}
  \item{corstr}{
a character string specifies the working correlation structure. The following are permitted: "independence", "exchangeable","ar1".
}
  \item{family}{
a description of the error distribution and link function to be used in the model. This is a character string naming a family function. The following are permitted: "gaussian", "binary", "poisson".
}
  \item{lag_level}{
specify how many lags of response y will affect missing probability.
}
}
\details{
Generate normal, poisson and binary longitudinal data based on the specified variance-covariance matrix.
}
\value{
\item{data}{simulated dataset}
\item{prob_miss}{missing percentage of the response y}
}

\author{
Cong Xu, Zheng Li and Ming Wang
}

\examples{
n=500
id=rep(1:n,each=3)
rho=1
phi=1
x=cbind(1,rnorm(length(id)))
beta=c(1,1)
x_mis=cbind(1,rnorm(length(id)))
para=c(1,1,1)
data_sim(id,rho,phi,x,beta,x_mis,para,"independent","gaussian",lag_level = 1)
}

\keyword{gee}
\keyword{simulate}% __ONLY ONE__ keyword per line
