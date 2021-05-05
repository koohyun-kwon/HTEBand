#' Lipschitz class worst-case bias
#'
#' Calculates the worst-case supremum bias for regression function value estimator under Lipschitz class
#'
#' @inheritParams K_fun
#' @param M bound on the first derivative
#'
#' @return a scalar worst-case bias value
#' @export
#'
bias_Lip <- function(x, t, M, kern, h){

  nmrt <- M * sum(K_fun(x, t, h, kern) * abs(x - t))
  dnmnt <- sum(K_fun(x, t, h, kern))

  return(nmrt / dnmnt)
}

#' Lipschitz class standard deviation
#'
#' Calculates the standard deviation for regression function value estimator under Lipschitz class
#'
#' @inheritParams K_fun
#' @inheritParams eps_hat
#'
#' @return a scalar standard deviation value
#' @export
var_Lip <- function(y, x, t, kern, h, deg, loo){

  sd.hat <- eps_hat(y, x, deg, kern, loo)

  nmrt <- sum(K_fun(x, t, h, kern)^2 * sd.hat^2)
  dnmnt <- sum(K_fun(x, t, h, kern))^2

  return(nmrt / dnmnt)
}

#' Lipschitz class true standard deviation
#'
#' Calculates the true standard deviation for regression function value estimator under Lipschitz class
#'
#' @inheritParams var_Lip
#' @param sd.true a vector of true conditional standard deviation values
#'
#' @return a scalar standard deviation value
#' @export
var_Lip_true <- function(x, t, kern, h, sd.true){

  nmrt <- sum(K_fun(x, t, h, kern)^2 * sd.true^2)
  dnmnt <- sum(K_fun(x, t, h, kern))^2

  return(nmrt / dnmnt)
}

# bw_Lip <- function(y, x, t, TE = FALSE, d = NULL, M, kern, alpha, se.initial, bw.eq = TRUE){
#
#   if(TE){
#
#     y.1 <- y[d == 1]
#     y.0 <- y[d == 0]
#     x.1 <- x[d == 1]
#     x.0 <- x[d == 0]
#
#     obj <- function(h.1, h.0){
#
#       bias <- M * (bias_Lip(x.1, t, M, kern, h.1) + bias_Lip(x.0, t, M, kern, h.0))
#       sd <-
#     }
#   }
#
#
# }
