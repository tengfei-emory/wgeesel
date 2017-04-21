\name{imps}
\alias{imps}
\non_function{}
\title{Inpatient Multidimensional Psychiatric Scale (IMPS)}
\usage{data(imps)}
\description{
The \code{imps} data frame has 1544 rows and 8 columns. The data is from National Institute of the Mental Health Schizophrenia Collaborative Study, where the effect of chlorpromazine, fluphenazine, or thioridazine treatment on the overall severity of the schizophrenia disorder is of interest.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{ID}{patient ID}
    \item{IMPS79}{the severity of the schizophrenia disorder (ranges from 0 to 7)}
    \item{Week}{the fixed visit time}
    \item{Drug}{the indicator of treatment (1: chlorphromazine, fluphenazine, or thioridazine treatment; 0: placebo)}
    \item{Sex}{the indicator of sex of the patients (1: male; 0: female) }    
    \item{R}{an indicator of the missingness (1: observed; 0: missing) }
    \item{Time}{square root of the Week covariate}
    \item{Y}{an indicator if IMPS \eqn{>=4} }
  }
}
\source{
Gibbons, R.D. and Hedeker, D., 1994. Application of random-effects probit regression models. \emph{Journal of consulting and clinical psychology}, 62(2), p.285.
}
\references{
Shen, C. W., & Chen, Y. H. (2012). Model selection for generalized estimating equations accommodating dropout missingness. \emph{Biometrics}, 68(4), 1046-1054.
}
\keyword{datasets}
