#' Residual calculation
#'
#' Calculates residuals after nonparametric regression using \code{\link[nprobust]{lprobust}}.
#'
#'
#' @param y vector of dependent variables
#' @param x vector of regressors
#' @param deg degree of local polynomial regression to be used
#' @param kern kernel used to calculate conditional variance function;
#' supports \code{"tri"}, \code{"epa"},
#' \code{"uni"}, and \code{"gau"}. Default is \code{kern = "epa"}.
#'
#' @return vector of residuals with the same length as \code{y}
#' @export
#'
#' @examples
#' x <- seq(from = -1, to = 1, length.out = 500)
#' y <- x^2 + stats::rnorm(500, 0, 0.1)
#' eps_hat(y, x, 1)
eps_hat <- function(y, x, deg, kern = "epa"){

  lp.res <- nprobust::lprobust(y, x, x, p = deg, kernel = kern)
  yhat.lp <- lp.res$Estimate[, "tau.us"]
  eps.hat <- y - yhat.lp

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
resid_calc <- function(y, x, d = NULL, deg, kern){

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

    resid.1[, j] <- eps_hat(y.1[, j], x.1, deg, kern)
    if(is.TE) resid.0[, j] <- eps_hat(y.0[, j], x.0, deg, kern)
  }

  return(list(resid.1 = resid.1, resid.0 = resid.0))
}
