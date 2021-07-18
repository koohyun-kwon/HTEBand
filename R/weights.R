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
#' @inheritParams w_get_Lip
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
w_get_Hol <- function(y, x, eval, C, level, kern = "triangular", se.initial = "EHW", se.method = "nn", J = 3,
                      TE = FALSE, d = NULL, bw.eq = TRUE){

  m <- length(eval)

  if(TE){

    if(kern == "triangular"){
      kern.rdh <- kern
      kern <- "tri"  # Naming convention is different
    }

    y.1 <- y[d == 1]
    y.0 <- y[d == 0]
    x.1 <- x[d == 1]
    x.0 <- x[d == 0]

    n.1 <- length(x.1)
    n.0 <- length(x.0)

    w.mat.1 <- matrix(0, n.1, m)
    w.mat.0 <- matrix(0, n.0, m)

    for(i in 1:m){

      h.opt <- bw_Lip(y, x, eval[i], TE, d, C, kern, 1 - level, bw.eq, p = 2)$h.opt
      d.1 <- RDHonest::LPPData(as.data.frame(cbind(y.1, x.1)), point = eval[i])
      d.0 <- RDHonest::LPPData(as.data.frame(cbind(y.0, x.0)), point = eval[i])

      w.mat.1[, i] <- RDHonest::NPRreg.fit(d.1, h.opt[1], kern = kern.rdh, se.method = se.method, J = J)$w
      w.mat.0[, i] <- RDHonest::NPRreg.fit(d.0, h.opt[2], kern = kern.rdh, se.method = se.method, J = J)$w
    }

    res <- list(w.mat.1 = w.mat.1, w.mat.0 = w.mat.0)

  }else{

    n <- length(y)
    w.mat <- matrix(0, nrow = n, ncol = m)

    for(i in 1:m){

      d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point = eval[i])

      bw.res.i <- RDHonest::NPROptBW.fit(d, M = C, kern = kern, opt.criterion = "OCI", alpha = 1 - level,
                                         beta = 0.5, se.initial = se.initial)

      w.res.i <- RDHonest::NPRreg.fit(d, bw.res.i$h[1], kern = kern, se.method = se.method, J = J)
      w.mat[, i] <- w.res.i$w
    }

    res <- w.mat
  }

  return(res)
}


#' Optimal weights for Lipschitz space
#'
#' Computes optimal weights and bandwidths for inference for nonparametric regression function values
#' under Lipschitz space.
#'
#' @inheritParams ci_reg_Lip
#' @param eval evaluation points
#'
#' @return a matrix of weights with each column corresponding to weights for each evaluation points,
#' or if \code{TE = TRUE}, a list of two such matrices, with \code{w.mat.1} corresponding to that of the
#' treated group and \code{w.mat.0} corresponding to that of the control group.
#' @export
#'
#' @examples
#' x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' eval <- seq(from = -0.9, to = 0.9, length.out = 5)
#' w_get_Lip(y, x, eval, 1, 0.95)
w_get_Lip <- function(y, x, eval, C, level, TE = FALSE, d = NULL, kern = "tri", bw.eq = TRUE){

  m <- length(eval)

  if(TE == TRUE){

    x.1 <- x[d == 1]
    x.0 <- x[d == 0]

    n.1 <- length(x.1)
    n.0 <- length(x.0)

    w.mat.1 <- matrix(0, n.1, m)
    w.mat.0 <- matrix(0, n.0, m)

    for(i in 1:m){

      h.opt <- bw_Lip(y, x, eval[i], TE, d, C, kern, 1 - level, bw.eq)$h.opt
      w.mat.1[, i] <- K_fun(x.1, eval[i], h.opt[1], kern) / sum(K_fun(x.1, eval[i], h.opt[1], kern))
      w.mat.0[, i] <- K_fun(x.0, eval[i], h.opt[2], kern) / sum(K_fun(x.0, eval[i], h.opt[2], kern))
    }

    res <- list(w.mat.1 = w.mat.1, w.mat.0 = w.mat.0)

  }else{

    n <- length(y)
    w.mat <- matrix(0, nrow = n, ncol = m)

    for(i in 1:m){

      h.opt <- bw_Lip(y, x, eval[i], TE, d, C, kern, 1 - level, bw.eq)$h.opt
      w.mat[, i] <- K_fun(x, eval[i], h.opt, kern) / sum(K_fun(x, eval[i], h.opt, kern))
    }

    res <- w.mat
  }

  return(res)
}
