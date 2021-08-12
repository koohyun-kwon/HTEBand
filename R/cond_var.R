#' Residual calculation
#'
#' Calculates residuals after nonparametric regression using \code{\link[nprobust]{lprobust}}
#' or \code{\link[locpol]{locPolSmootherC}}.
#'
#' For now, only \code{kern = "tri"} and \code{kern = "epa"} are supported when
#' \code{var.reg = "locpol"}.
#'
#'
#' @param y vector of dependent variables
#' @param x vector of regressors
#' @param deg degree of local polynomial regression to be used
#' @param kern kernel used to calculate conditional variance function;
#' supports \code{"tri"}, \code{"epa"},
#' \code{"uni"}, and \code{"gau"}. Default is \code{kern = "epa"}.
#' @param var.reg nonparametric regression method used to calculate residuals;
#' either \code{"npr"} (standing for \code{nprobust} package) or \code{"locpol"}
#' (standing for \code{locpol} package). Default is \code{var.reg = "npr"}.
#' @param n.max maximum number of observations allowed where the regression function estimator
#' is evaluated when \code{var.reg = "npr"}.
#' If \code{length(x) > n.max}, the regression function estimator
#' is evaluated over a grid of the length \code{n.max}, and the regression function estimate over
#' \code{x} is obtained by linear interpolation. The default is \code{n.max = 500}.
#'
#' @return vector of residuals with the same length as \code{y}
#' @export
#'
#' @examples
#' x <- seq(from = -1, to = 1, length.out = 500)
#' y <- x^2 + stats::rnorm(500, 0, 0.1)
#' eps_hat(y, x, 1)
#' eps_hat(y, x, 1, var.reg = "locpol")
eps_hat <- function(y, x, deg, kern = "epa", var.reg = "npr", n.max = 500){

  if(var.reg == "npr"){
    if(length(x) > n.max){
      x.eval <- seq(from = min(x), to = max(x), length.out = n.max)
      lp.res <- nprobust::lprobust(y, x, x.eval, p = deg, kernel = kern, vce = "hc0")
      yhat.grid <- lp.res$Estimate[, "tau.us"]
      yhat.lp <- stats::approx(x.eval, yhat.grid, x)$y
    }else{
      lp.res <- nprobust::lprobust(y, x, x, p = deg, kernel = kern)
      yhat.lp <- lp.res$Estimate[, "tau.us"]
    }
    eps.hat <- y - yhat.lp

  }else if(var.reg == "locpol"){

    if(kern == "tri"){
      K <- locpol::TrianK
    }else if(kern == "epa"){
      K <- locpol::EpaK
    }

    bw <- locpol::thumbBw(x, y, deg, K)
    locpol.res <- locpol::locPolSmootherC(x, y, x, bw, deg, K)
    y.hat <- locpol.res$beta0
    eps.hat <- y - y.hat
  }

  return(eps.hat)
}


#' Residual calculation function
#'
#' Calculates residuals corresponding to each observations.
#'
#' @inheritParams opt_w
#'
#' @return a list of residuals corresponding to treated individuals (\code{resid.1}) and control individuals (\code{resid.0});
#' if there is no control group, \code{resid.0 = 0} is returned.
#' @export
resid_calc <- function(y, x, d = NULL, deg, kern = "epa", var.reg = "npr"){

  is.TE <- sum(d == 0) > 0
  if(!is.TE) d <- rep(1, length(y))

  y <- v_to_m(y)
  k <- ncol(y)

  y.1 <- v_to_m(y[d == 1, ])
  y.0 <-
    if(sum(d == 0) > 0){
      v_to_m(y[d == 0, ])
    }else{
      v_to_m(0)
    }

  x.1 <- x[d == 1]
  x.0 <-
    if(is.TE){
      x[d == 0]
    }else{
      0
    }

  n.1 <- length(x.1)
  n.0 <- length(x.0)

  resid.1 <- matrix(0, nrow = n.1, ncol = k)
  resid.0 <- matrix(0, nrow = n.0, ncol = k)

  for(j in 1:k){

    resid.1[, j] <- eps_hat(y.1[, j], x.1, deg, kern, var.reg)
    if(is.TE) resid.0[, j] <- eps_hat(y.0[, j], x.0, deg, kern, var.reg)
  }

  return(list(resid.1 = resid.1, resid.0 = resid.0))
}
