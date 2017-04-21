\name{ylag}
\alias{ylag}
\encoding{utf-8}
\title{
Generate subject-level lagged response
}
\description{
The function generates lagged response for each subject.
}
\usage{
ylag(id,y,lag)
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

}
\value{
return a vector of lagged y
}

\author{
Cong Xu, Zheng Li and Ming Wang
}

\examples{
id=rep(c(1:20),each=3)
lag=1
y=rnorm(length(id))
ylag(id,y,1)
}

\keyword{gee}
\keyword{lag}% __ONLY ONE__ keyword per line
