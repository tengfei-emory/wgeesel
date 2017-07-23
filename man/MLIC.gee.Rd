\name{MLIC.gee}
\alias{MLIC.gee}
\encoding{utf-8}
\title{
MLIC and MLICC for Weighted GEE
}
\description{
Calculate the MLIC (missing longitudinal information
criterion) for selection of mean model, and the MLICC (missing longitudinal information correlation
criterion) for selection of working correlation structure, based on the expected quadratic loss and the WGEE.
}
\usage{
MLIC.gee(object,object_full)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted model object of class \code{"wgee"}.
}
 \item{object_full}{
a fitted model object of class \code{"wgee"}: the largest candidate model under consideration to be fitted. 
}

}
 
\value{
Return a data frame of MLIC, MLICC and Wquad_loss.
}

\author{
Cong Xu, Zheng Li and Ming Wang
}

\references{
Robins, J.M., Rotnitzky, A. and Zhao, L.P., 1995. Analysis of semiparametric regression models for repeated outcomes in the presence of missing data. \emph{Journal of the American Statistical Association}, 90(429), pp.106-121.

Shen, C.W. and Chen, Y.H., 2012. Model selection for generalized estimating equations accommodating dropout missingness. \emph{Biometrics}, 68(4), pp.1046-1054.

Shen, C.W. and Chen, Y.H., 2013. Model selection of generalized estimating equations with multiply imputed longitudinal data. \emph{Biometrical Journal}, 55(6), pp.899-911.
}
\note{
MLIC and MLICC model selection criterion for longitudinal data criterion with dropouts or monotone missingness under the assumption of MAR.
}

\seealso{
\code{\link{wgee}}
}

\examples{
data(imps)
fit1 <- wgee (Y ~ Drug+Sex+Time,data=imps,id=imps$ID,family="binomial",
             corstr="exchangeable",scale=NULL,mismodel= R ~ Drug+Time)

fit_f <- wgee (Y ~ Drug+Sex+Time+Time*Sex+Time*Drug,data=imps,id=imps$ID, family="binomial",
              corstr="exchangeable",scale=NULL,mismodel= R ~ Drug+Time)
###not run#####
##MLIC.gee(fit1,fit_f)

}

\keyword{MLIC}
\keyword{wgee}% __ONLY ONE__ keyword per line
