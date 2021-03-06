% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_const.R
\name{ci_reg_Hol}
\alias{ci_reg_Hol}
\title{Confidence interval for regression function value for Hölder space}
\usage{
ci_reg_Hol(
  y,
  x,
  point,
  C,
  level,
  kern = "triangular",
  se.initial = "EHW",
  se.method = "nn",
  J = 3,
  cv = NULL,
  TE = FALSE,
  d = NULL,
  bw.eq = TRUE
)
}
\arguments{
\item{y}{a vector of dependent variable}

\item{x}{a vector of regressor}

\item{point}{point where the regression function value would be evaluated}

\item{C}{bound on the second derivative}

\item{level}{confidence level of each one-sided confidence intervals}

\item{kern}{specifies kernel function used in the local regression; default = \code{"triangular"}.
See \code{\link[RDHonest]{NPROptBW.fit}} in \code{RDHonest} package for a list of kernels available.}

\item{se.initial}{method for estimating initial variance for computing optimal bandwidt; default = \code{"EHW"}.
See \code{\link[RDHonest]{NPROptBW.fit}} in \code{RDHonest} package for a list of method available.}

\item{se.method}{methods for estimating standard error of estimate; default = \code{"nn"}.
See \code{\link[RDHonest]{NPRreg.fit}} in \code{RDHonest} package for a list of method available.}

\item{J}{number of nearest neighbors, if "nn" is specified in se.method.}

\item{cv}{supplied value of critical value to be used in constructing confidence interval;
default is \code{cv = NULL}.}

\item{TE}{logical specifying whether there are treatment and control groups.}

\item{d}{a vector of indicator variables specifying treatment and control group status;
relevant only when \code{TE = TRUE}.}

\item{bw.eq}{if \code{TRUE}, the same bandwidths are used for estimators for treatment and control groups;
relevant only when \code{TE = TRUE}.}
}
\value{
a vector of lower and upper ends of the confidence interval  and a pair of bandwidths used for
the treatment and control groups.
}
\description{
Constructs a confidence interval for regression function value for Hölder space based on
the optimal one-sided procedure of Armstrong and Kolesár (2020).
}
\examples{
x <- stats::runif(500, min = -1, max = 1)
y <- x + rnorm(500, 0, 1/4)
ci_reg_Hol(y, x, 1/2, 1, 0.99)
}
\references{
Armstrong, Timothy B., and Michal Kolesár. 2020.
"Simple and Honest Confidence Intervals in Nonparametric Regression." Quantitative Economics 11 (1): 1–39.
}
