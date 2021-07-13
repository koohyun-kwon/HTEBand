#' Residual calculation
#'
#' Calculates residuals after nonparametric regression using \code{\link[locpol]{locPolSmootherC}}.
#'
#' The bandwidth is calculated using the rule-of-thumb method of Fan and Gijbels (1996)
#' implemented by \code{\link[locpol]{thumbBw}}.
#'
#' @param y vector of dependent variables
#' @param x vector of regressors
#' @param deg degree of local polynomial regression to be used
#' @param kern name of the kernel. Supports \code{"triangle"} and \code{"epa"}.
#' @param loo If \code{TRUE}, the leave-one-out version is used; see \code{\link[locpol]{locPolSmootherC}}.
#'
#' @return vector of residuals with the same length as \code{y}
#' @export
#'
#' @examples
#' x <- seq(from = -1, to = 1, length.out = 500)
#' y <- x^2 + stats::rnorm(500, 0, 0.1)
#' eps_hat(y, x, 0, "triangle", FALSE)
#' eps_hat(y, x, 1, "epa", TRUE)
eps_hat <- function(y, x, deg, kern, loo){


  if(kern == "triangle"){

    K <- locpol::TrianK
  }else if(kern == "epa"){

    K <- locpol::EpaK
  }

  bw <- locpol::thumbBw(x, y, deg, K)

  if(loo){
    locpol.res <- locpol::looLocPolSmootherC(x, y, bw, deg, K)
  }else{
    locpol.res <- locpol::locPolSmootherC(x, y, x, bw, deg, K)
  }

  y.hat <- locpol.res$beta0
  eps.hat <- y - y.hat

  return(eps.hat)
}
