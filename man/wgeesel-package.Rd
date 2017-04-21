\name{wgeesel-package}
\alias{wgeesel-package}
\alias{wgeesel}
\docType{package}
\encoding{utf-8}
\title{
Weighted Generalized Estimating Equations and Model Selection
}
\description{
Weighted Generalized estimating equations (WGEE) is an extension of generalized linear models to longitudinal or clustered data by incorporating the correlation within-cluster when data is missing at random (MAR). The parameters in mean, scale, correlation structures are estimate based on quasi-likelihood. The package \pkg{wgeesel} also contains model selection criteria for variable selection in the mean model and for the selection of a working correlation structure in longitudinal data with dropout or monotone missingness using WGEE.
}

\details{

The collection of functions includes:
\describe{
\item{\code{\link{wgee}}}{estimates parameters based on WGEE in
mean, scale, and correlation structures, through mean link, scale link, and correlation link}.

\item{\code{\link{QIC.gee}}, \code{\link{MQIC.gee}}, \code{\link{RJ.gee}}}{calculate the QIC (QIC\eqn{_{u}}{u}), MQIC (MQIC\eqn{_{u}}{u}), Rotnitzky-Jewell criteria for variable selection in the mean model and/or selection of a working correlation structure in GEE (unbalanced data is allowed).}

\item{\code{\link{MLIC.gee}}, \code{\link{QICW.gee}}}{calculate the MLIC (MLICC) and QICW\eqn{_{r}}{r} (QICW\eqn{_{p}}{p}) for variable selection in the mean model and the selection of a working correlation structure in WGEE, which can accommodate dropout missing at random (MAR).} 

\item{\code{\link{data_sim}}}{can simulate longitudinal response data in different distribution (gaussian, binomial, poisson) with drop missingness.}
}

For a complete list of functions, use \code{library(help = "wgeesel")}.
}

\author{
Cong Xu \email{congxu17@gmail.com}, 
Zheng Li \email{zheng.li@outlook.com}, 
Ming Wang \email{mwang@phs.psu.edu}


Maintainer:
Zheng Li \email{zheng.li@outlook.com}
}

\references{
Liang, K.Y. and Zeger, S.L., 1986. Longitudinal data analysis using generalized linear models. \emph{Biometrika}, pp.13-22.

Preisser, J.S., Lohman, K.K. and Rathouz, P.J., 2002. Performance of weighted estimating equations for longitudinal binary data with drop-outs missing at random. \emph{Statistics in medicine}, 21(20), pp.3035-3054.

Robins, J.M., Rotnitzky, A. and Zhao, L.P., 1995. Analysis of semiparametric regression models for repeated outcomes in the presence of missing data. \emph{Journal of the American Statistical Association}, 90(429), pp.106-121.

Shen, C. W., & Chen, Y. H. (2012). Model selection for generalized estimating equations accommodating dropout missingness. \emph{Biometrics}, 68(4), 1046-1054.

Wang, M., 2014. Generalized Estimating Equations in Longitudinal Data Analysis: A Review and Recent Developments. \emph{Advances in Statistics}, 2014.
}

\keyword{gee}
\keyword{wgee}
\keyword{package}
\keyword{qic}

\seealso{
GEE methods exist for \code{\link[geepack]{geeglm}} (\pkg{geepack})
}

\examples{
data(imps)

fit <- wgee (Y ~ Drug+Sex+Time,data=imps,id=imps$ID,family="binomial",
  corstr="exchangeable",scale=NULL,mismodel= R ~ Drug+Time)

#MLIC.gee(Y ~ Drug+Sex+Time, model_full=Y ~ Drug+Sex+Time+Time*Sex+Time*Drug,
#  mismodel= R ~ Drug+Time, data=imps,id=imps$ID, family="binomial",
#  corstr="exchangeable")
  
QICW.gee(Y ~ Drug+Sex+Time, mismodel= R ~ Drug+Time, data=imps, id=imps$ID,
  family="binomial", corstr="exchangeable")
}
