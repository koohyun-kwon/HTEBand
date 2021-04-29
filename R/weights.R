#' Optimal weights for Hölder space
#'
#' Computes optimal weights and bandwidths for inference for nonparametric regression function values
#' under Hölder space, with an aid of \code{RDHonest} package.
#'
#' @param y vector of dependent variables
#' @param x vector of regressors
#' @param eval evaluation points
#' @param C bound on the second derivative
#' @param level confidence level
#' @param kern specifies kernel function used in the local regression; default = \code{"triangular"}.
#' See \code{\link[RDHonest]{NPROptBW.fit}} in \code{RDHonest} package for a list of kernels available.
#' @param se.initial method for estimating initial variance for computing optimal bandwidt; default = \code{"EHW"}.
#' See \code{\link[RDHonest]{NPROptBW.fit}} in \code{RDHonest} package for a list of method available.
#' @param se.method methods for estimating standard error of estimate; default = \code{"nn"}.
#' See \code{\link[RDHonest]{NPRreg.fit}} in \code{RDHonest} package for a list of method available.
#' @param J number of nearest neighbors, if "nn" is specified in se.method.
#'
#' @return list with components
#'
#' \describe{
#'     \item{w.mat}{matrix of optimal weigths with a dimension \code{length(y)} by \code{length(eval)}
#'     corresponding to the evaluation points}
#'     \item{bw.vec}{vector of optimal bandwidths corresponding to the evaluation points}
#'
#' }
#' @export
#'
#' @examples x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' eval <- seq(from = -0.9, to = 0.9, length.out = 5)
#' w_get_Hol(y, x, eval, 1, 0.95)
w_get_Hol <- function(y, x, eval, C, level, kern = "triangular", se.initial = "EHW", se.method = "nn", J = 3){

  n <- length(y)
  m <- length(eval)
  bw.vec <- numeric(m)
  w.mat <- matrix(0, nrow = n, ncol = m)

  for(i in 1:m){

    d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point = eval[i])

    bw.res.i <- RDHonest::NPROptBW.fit(d, M = C, kern = kern, opt.criterion = "OCI", alpha = 1 - level,
                                      beta = 0.5, se.initial = se.initial)
    bw.vec[i] <- bw.res.i$h[1]

    w.res.i <- RDHonest::NPRreg.fit(d, bw.vec[i], kern = kern, se.method = se.method, J = J)
    w.mat[, i] <- w.res.i$w
  }

  res <- list(w.mat = w.mat, bw.vec = bw.vec)

  return(res)
}
