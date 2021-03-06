% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_const.R
\name{ci_reg_Lip}
\alias{ci_reg_Lip}
\title{Confidence interval for regression function value for Lipschitz space}
\usage{
ci_reg_Lip(
  y,
  x,
  point,
  C,
  level,
  TE = FALSE,
  d = NULL,
  kern = "tri",
  bw.eq = TRUE,
  se.method = "resid",
  cv = NULL
)
}
\arguments{
\item{y}{vector of dependent variables}

\item{x}{a vector of regressors.}

\item{point}{point where the regression function value would be evaluated}

\item{C}{bound on the second derivative}

\item{level}{confidence level of each one-sided confidence intervals}

\item{TE}{logical specifying whether there are treatment and control groups.}

\item{d}{a vector of indicator variables specifying treatment and control group status;
relevant only when \code{TE = TRUE}.}

\item{kern}{a string for kernel name; currently \code{"tri"} is supported.}

\item{bw.eq}{if \code{TRUE}, the same bandwidths are used for estimators for treatment and control groups;
relevant only when \code{TE = TRUE}.}

\item{se.method}{methods for estimating standard error of estimate; currently,
only "resid" is supported.}

\item{cv}{supplied value of critical value to be used in constructing confidence interval;
default is \code{cv = NULL}.}
}
\value{
a vector of lower and upper ends of confidence interval, and a pair of bandwidths used for
the treatment and control groups; if \code{TE = FALSE}, those two bandwidths have the same value.
}
\description{
Constructs a confidence interval for regression function value for Lipschitz space;
it can be used to deal with both regression function value itself and
difference between regression function values of treatment and control groups.
}
\examples{
x <- stats::runif(500, min = -1, max = 1)
y <- x + rnorm(500, 0, 1/4)
ci_reg_Lip(y, x, 1/2, 1, 0.99)
}
