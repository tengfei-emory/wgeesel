\name{ylag}
\alias{ylag}
\encoding{utf-8}
\title{
Generate subject-level lagged response
}
\description{
Suppose you have longitudinal response y and it's subject id. This function generates lagged y for each subject.
}
\usage{
ylag(id,y,lag,na=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
subject id
}
  \item{y}{
response
}

  \item{lag}{
how many lags for y
}
  \item{na}{
logical; whether remain NAs in the response when lag \eqn{>1}. Default to FALSE and output 0s instead of NAs.
}
}
\value{
return a vector of lagged y
}

\author{
Cong Xu, Zheng Li and Ming Wang
}

\examples{
id <- rep(c(1:20),each=3)
y <- rnorm(length(id))
ylag(id,y,1) #lag=1
ylag(id,y,2,na=FALSE) #lag=2
}

\keyword{gee}
\keyword{lag}% __ONLY ONE__ keyword per line
