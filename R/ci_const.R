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
#' @return a vector of lower and upper ends of the confidence interval
#' @export
#' @references Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression." Quantitative Economics 11 (1): 1–39.
#'
#' @examples
#' x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' ci_reg_Hol(y, x, 1/2, 1, 0.99)
ci_reg_Hol <- function(y, x, point, C, level, kern = "triangular", se.initial = "EHW", se.method = "nn", J = 3){

  d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point)

  ci.res <- RDHonest::NPRHonest.fit(d = d, M = C, kern = kern, opt.criterion = "OCI",
                                    alpha = 1 - level, beta = 0.5, se.method = se.method, J = J,
                                    se.initial = se.initial)

  return(c(ci.res$lower, ci.res$upper))
}
