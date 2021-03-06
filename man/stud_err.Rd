% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile.R
\name{stud_err}
\alias{stud_err}
\title{Studentized error variable}
\usage{
stud_err(w.1, w.0, resid.1, resid.0, omega.1, omega.0, T.grad)
}
\arguments{
\item{w.1}{A \code{n_1} by \code{k} dimensional matrix of weight values corresponding to treated observations}

\item{w.0}{A \code{n_0} by \code{k} dimensional matrix of weight values corresponding to control observations}

\item{resid.1}{A \code{n_1} by \code{k} dimensional matrix of residuals corresponding to treated observations}

\item{resid.0}{A \code{n_0} by \code{k} dimensional matrix of residuals corresponding to control observations}

\item{omega.1}{A \code{n_1} by \code{k} by \code{k} array of conditional variance corresponding to treated observations}

\item{omega.0}{A \code{n_0} by \code{k} by \code{k} array of conditional variance corresponding to control observations}

\item{T.grad}{A \code{k} dimensional gradient vector of \code{T_t f}}
}
\value{
A studentized scalar estimation error term
}
\description{
Calculates the studentized error term given weights, residuals, and conditional variance matrix.
}
\examples{
stud_err(rep(1/50, 50), rep(1/50, 50), stats::rnorm(50),
stats::rnorm(50), rep(1, 50), rep(1, 50), 1)
w.1 <- w.0 <- matrix(rep(1/50, 100), ncol = 2)
resid.1 <- resid.0 <- matrix(stats::rnorm(100), ncol = 2)
omega.1 <- omega.0 <- array(rep(c(1,0,0,1), each = 50), dim = c(50, 2, 2))
T.grad <- c(1, 1)
stud_err(w.1, w.0, resid.1, resid.0, omega.1, omega.0, T.grad)
}
