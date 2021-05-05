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

#' Lipschitz class variance
#'
#' Calculates the variance for regression function value estimator under Lipschitz class
#'
#' @inheritParams K_fun
#' @inheritParams eps_hat
#'
#' @return a scalar variance value
#' @export
var_Lip <- function(y, x, t, kern, h, deg, loo){

  sd.hat <- eps_hat(y, x, deg, kern, loo)

  nmrt <- sum(K_fun(x, t, h, kern)^2 * sd.hat^2)
  dnmnt <- sum(K_fun(x, t, h, kern))^2

  return(nmrt / dnmnt)
}

#' Lipschitz class true variance
#'
#' Calculates the true variance for regression function value estimator under Lipschitz class
#'
#' @inheritParams var_Lip
#' @param sd.true a vector of true conditional standard deviation values
#'
#' @return a scalar variance value
#' @export
var_Lip_true <- function(x, t, kern, h, sd.true){

  nmrt <- sum(K_fun(x, t, h, kern)^2 * sd.true^2)
  dnmnt <- sum(K_fun(x, t, h, kern))^2

  return(nmrt / dnmnt)
}

#' Optimal bandwidths under Lipschitz class
#'
#' Calculates the optimal bandwidths for regression function value estimator under Lipschitz class
#'
#' @inheritParams bias_Lip
#' @inheritParams var_Lip
#' @param TE logical specifying whether there are treatment and control groups
#' @param d a vector of indicator variables specifying treatment and control group status;
#' relevant only when \code{TE = TRUE}.
#' @param alpha determines confidence level \code{1 - alpha}
#' @param bw.eq if \code{TRUE}, the same bandwidths are used for estimators for treatment and control groups;
#' relevant only when \code{TE = TRUE}.
#'
#' @return a list with the following components
#' \describe{
#'
#' \item{h.opt}{the optimal bandwidth; when \code{TE = TRUE}, return two bandwidths
#' for estimators for treatment and control groups.}
#'
#' \item{hl.opt}{the optimal half-length when the optimal bandwidth(s) is used.}
#' }
#' @export
bw_Lip <- function(y, x, t, TE = FALSE, d = NULL, M, kern, alpha, bw.eq = TRUE,
                   deg, loo){

  if(TE == TRUE){

    y.1 <- y[d == 1]
    y.0 <- y[d == 0]
    x.1 <- x[d == 1]
    x.0 <- x[d == 0]

    obj <- function(h){

      h.1 <- abs(h[1]) # optim() might evaluate negative bandwidths
      h.0 <- abs(h[2])
      bias <- M * (bias_Lip(x.1, t, M, kern, h.1) + bias_Lip(x.0, t, M, kern, h.0))
      sd <- sqrt(var_Lip(y.1, x.1, t, kern, h.1, deg, loo) + var_Lip(y.0, x.0, t, kern, h.0, deg, loo))
      c <- stats::qnorm(1 - alpha) / 2
      return(bias + c * sd)
    }

    if(bw.eq == FALSE){

      h.1.init <- max(abs(x.1 - t)) / 2
      h.0.init <- max(abs(x.0 - t)) / 2

      opt.res <- stats::optim(c(h.1.init, h.0.init), obj)
      h.opt <- abs(opt.res$par) # optim() might evaluate negative bandwidths
      hl.opt <- opt.res$value  # half-length
    }else{

      obj.eq <- function(h){
        obj(c(h, h))
      }

      h.min <- max(sort(unique(abs(x.1 - t)))[2], sort(unique(abs(x.0 - t)))[2])
      h.max <- max(abs(x.1 - t), abs(x.0 - t))

      opt.res <- stats::optimize(obj.eq, c(h.min, h.max), tol = .Machine$double.eps^0.75)
      h.opt <- rep(opt.res$minimum, 2)
      hl.opt <- opt.res$objective
    }

  }else{

    obj.1 <- function(h){

      bias <- M * bias_Lip(x, t, M, kern, h)
      sd <- sqrt(var_Lip(y, x, t, kern, h, deg, loo))
      c <- stats::qnorm(1 - alpha) / 2
      return(bias + c * sd)
    }

    h.min <- sort(unique(abs(x - t)))[2]
    h.max <- abs(x - t)

    opt.res <- stats::optimize(obj.1, c(h.min, h.max), tol = .Machine$double.eps^0.75)
    h.opt <- opt.res$minimum
    hl.opt <- opt.res$objective
  }

  res <- list(h.opt = h.opt, hl.opt = hl.opt)
  return(res)
}
