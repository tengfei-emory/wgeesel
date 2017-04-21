\name{seizure}
\alias{seizure}
\non_function{}
\title{Epiliptic Seizures}
\usage{data(seizure)}
\description{
The dataset has the number of epiliptic seizures in each of four two-week 
intervals,and in a baseline eight-week inverval, for treatment and control 
groups with a total of 59 individuals. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y1}{the number of epiliptic seizures in the 1st 2-week
      interval}
    \item{y2}{the number of epiliptic seizures in the 2nd 2-week
      interval}
    \item{y3}{the number of epiliptic seizures in the 3rd 2-week
      interval}
    \item{y4}{the number of epiliptic seizures in the 4th 2-week
      interval}    
    \item{trt}{an indicator of treatment}
    \item{base}{the number of epilitic seizures in a baseline 8-week
      interval}
    \item{age}{a numeric vector of subject age}
  }
}
\source{
  Thall, P.F. and Vail S.C. (1990)
  Some covariance models for longitudinal count data with
  overdispersion. \emph{Biometrics} \bold{46}: 657--671.
}
\references{
  Diggle, P.J., Liang, K.Y., and Zeger, S.L. (1994)
  Analysis of Longitudinal Data. Clarendon Press.
}
\seealso{
\code{\link[geepack]{seizure}} (\pkg{geepack})
}
\keyword{datasets}
