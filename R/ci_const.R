#' Confidence interval for regression function value for Hölder space
#'
#' Constructs a confidence interval for regression function value for Hölder space based on
#' the optimal one-sided procedure of Armstrong and Kolesár (2020).
#'
#' @param y a vector of dependent variable
#' @param x a vector of regressor
#' @param point point where the regression function value would be evaluated
#' @param C bound on the second derivative
#' @param level confidence level of each one-sided confidence intervals
#' @inheritParams w_get_Hol
#'
#' @return a vector of lower and upper ends of the confidence interval  and a pair of bandwidths used for
#' the treatment and control groups.
#' @export
#' @references Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression." Quantitative Economics 11 (1): 1–39.
#'
#' @examples
#' x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' ci_reg_Hol(y, x, 1/2, 1, 0.99)
ci_reg_Hol <- function(y, x, point, C, level, kern = "triangular", se.initial = "EHW",
                       se.method = "nn", J = 3){

  d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point)

  ci.res <- RDHonest::NPRHonest.fit(d = d, M = C, kern = kern, opt.criterion = "OCI",
                                    alpha = 1 - level, beta = 0.5, se.method = se.method, J = J,
                                    se.initial = se.initial)

  maxbias <- ci.res$maxbias
  sd <- unname(ci.res$sd)
  c <- stats::qnorm(level)/2
  hl <- maxbias + sd * c
  ci.lower <- ci.res$estimate - hl
  ci.upper <- ci.res$estimate + hl

  return(c(ci.lower, ci.upper, ci.res$hp, ci.res$hp))
}

#' Confidence interval for regression function value for Lipschitz space
#'
#' Constructs a confidence interval for regression function value for Lipschitz space;
#' it can be used to deal with both regression function value itself and
#' difference between regression function values of treatment and control groups.
#'
#' @inheritParams bw_Lip
#' @inheritParams ci_reg_Hol
#' @param se.method methods for estimating standard error of estimate; currently,
#' only "resid" is supported.
#'
#' @return a vector of lower and upper ends of confidence interval, and a pair of bandwidths used for
#' the treatment and control groups; if \code{TE = FALSE}, those two bandwidths have the same value.
#' @export
#'
#' @examples
#' x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' ci_reg_Lip(y, x, 1/2, 1, 0.99)
ci_reg_Lip <- function(y, x, point, C, level, TE = FALSE, d = NULL, kern = "tri",
                       bw.eq = TRUE, deg = 0, se.method = "resid"){

  opt.res <- bw_Lip(y, x, point, TE, d, C, kern, 1 - level, bw.eq, deg)
  h.opt <- opt.res$h.opt

  if(TE == TRUE){

    y.1 <- y[d == 1]
    y.0 <- y[d == 0]
    x.1 <- x[d == 1]
    x.0 <- x[d == 0]

    est.1 <- sum(K_fun(x.1, point, h.opt[1], kern) * y.1) / sum(K_fun(x.1, point, h.opt[1], kern))
    est.0 <- sum(K_fun(x.0, point, h.opt[2], kern) * y.0) / sum(K_fun(x.0, point, h.opt[2], kern))
    est <- est.1 - est.0
  }else{

    est <- sum(K_fun(x, point, h.opt, kern) * y) / sum(K_fun(x, point, h.opt, kern))
    h.opt <- rep(h.opt, 2)  # To make the returned value consistent
  }

  if(se.method == "resid"){

    ci.lower <- est - opt.res$hl.opt
    ci.upper <- est + opt.res$hl.opt
  }

  return(c(ci.lower, ci.upper, h.opt))
}
