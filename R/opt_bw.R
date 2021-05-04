#' Lipschitz class worst-case bias
#'
#' Calculates the worst-case bias for regression function values under Lipschitz class
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

# bw_Lip <- function(y, x, t, TE, d, M, kern, alpha, beta, se.initial, bw.eq){
#
#
# }
