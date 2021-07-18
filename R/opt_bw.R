#' Lipschitz class worst-case bias
#'
#' Calculates the worst-case supremum bias for regression function value estimator under Lipschitz class.
#'
#' @inheritParams K_fun
#' @param M bound on the first derivative
#'
#' @return a scalar worst-case bias value
#' @export
#'
bias_Lip <- function(x, t, M, kern, h){

  h.min <- min(abs(x - t))

  if(h < 0){
    stop("Negative bandwidth is not allowed")
  }else if(h <= h.min){

    # As h approaches h.min from above, the bias approaches M * h.min
    # This is in order to preserve the asymptotic bias representation that bias = h * C
    # for some constant C > 0.
    res <- M * h
  }else{

    nmrt <- M * sum(K_fun(x, t, h, kern) * abs(x - t))
    dnmnt <- sum(K_fun(x, t, h, kern))

    res <- nmrt/dnmnt
  }

  return(res)
}


#' Hölder class worst-case bias
#'
#' Calculates the worst-case supremum bias for regression function value estimator
#' under Hölder class.
#'
#' @inheritParams bias_Lip
#' @param M bound on the second derivative
#'
#' @return a scalar worst-case bias value
#' @export
#'
bias_Hol <- function(x, t, M, kern, h){

  # Construct LPPData() object; we are taking y = x, which don't matter
  d <- RDHonest::LPPData(as.data.frame(cbind(x, x)), point = t)
  d <- RDHonest::NPRPrelimVar.fit(d, se.initial = "EHW") # Add conditional variance values

  if(kern == "tri"){
    kern.rdh = "triangular"  # Naming convention is different from mine
  }

  r1 <- RDHonest::NPRreg.fit(d, h, kern.rdh, order = 1, se.method = "supplied.var",
                             TRUE) # NPRreg.fit() accepts scalar h

  w <- r1$w  # local linear weight values
  wt <- w[w != 0]
  xx <- d$X[w != 0]
  nobs <- length(wt)

  if (nobs == 0) bias <- sqrt(.Machine$double.xmax/10)
  # nobs = 0 when h is very small;
  # by taking very large bias, prevents very small h from being chosen during optimization

  w2p <- function(s) abs(sum((wt * (xx - s))[xx >= s]))
  w2m <- function(s) abs(sum((wt * (s - xx))[xx <= s]))
  bp <- stats::integrate(function(s) vapply(s, w2p,
                                            numeric(1)), 0, h)$value
  bm <- stats::integrate(function(s) vapply(s, w2m,
                                            numeric(1)), -h, 0)$value
  bias <- M * (bp + bm)

  return(bias)
}

#' Lipschitz class variance
#'
#' Calculates the variance for regression function value estimator under Lipschitz class
#'
#' @inheritParams K_fun
#' @inheritParams eps_hat
#' @param sd.homo logical indicating whether the variance would be estimated under locally homoskedastic variance
#' assumption; the default is \code{sd.homo = TRUE}.
#'
#' @return a scalar variance value
#' @export
var_Lip <- function(y, x, t, kern, h, deg = NULL, sd.homo = TRUE){

  h.min <- min(abs(x - t))

  if(sd.homo == TRUE){

    d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point = t)
    d <- RDHonest::NPRPrelimVar.fit(d, se.initial = "EHW")
    sd.hat <- sqrt(d$sigma2)
  }else{

    sd.hat <- eps_hat(y, x, deg)  # Not used right now
  }

  if(h < 0){

    stop("Negative bandwidth is not allowed")

  }else if(h <= h.min){

    # As h approaches h.min from above, the variance approaches sigma^2_{i(min)}
    # This is in order to preserve the asymptotic bias representation that
    # var = C_n / h for some constant C_n > 0.

    min.ind <- which.min(abs(x - t))
    res <- h.min * sd.hat[min.ind]^2 / h
  }else{

    nmrt <- sum(K_fun(x, t, h, kern)^2 * sd.hat^2)
    dnmnt <- sum(K_fun(x, t, h, kern))^2

    res <- nmrt / dnmnt
  }

  return(res)
}


#' Hölder class variance
#'
#' Calculates the variance for regression function value estimator under Hölder class
#'
#' @inheritParams K_fun
#' @inheritParams eps_hat
#'
#' @return a scalar variance value
#' @export
var_Hol <- function(y, x, t, kern, h){

  if(kern == "tri"){
    kern.rdh = "triangular"  # Naming convention is different from mine
  }

  # Construct LPPData() object;
  d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point = t)
  d <- RDHonest::NPRPrelimVar.fit(d, se.initial = "EHW") # Add conditional variance values

  r1 <- RDHonest::NPRreg.fit(d, h, kern.rdh, order = 1, se.method = "supplied.var",
                             TRUE) # NPRreg.fit() accepts scalar h

  return(r1$se["supplied.var"]^2)
}




#' Lipschitz class variance using residuals (obsolete now)
#'
#' Calculates the variance for regression function value estimator under Lipschitz class using residuals
#'
#' When \code{resid} corresponds to the true conditional standard deviations, this function calculates
#' the true variance value.
#'
#' @inheritParams var_Lip
#' @param resid a vector of true conditional standard deviation values
#'
#' @return a scalar variance value
#' @export
var_Lip_resid <- function(x, t, kern, h, resid){

  if(h <= 0){

    res <- 0
  }else{

    nmrt <- sum(K_fun(x, t, h, kern)^2 * resid^2)
    dnmnt <- sum(K_fun(x, t, h, kern))^2
    res <- nmrt/ dnmnt
  }

  return(res)
}

#' Optimal bandwidths under Lipschitz class
#'
#' Calculates the optimal bandwidths for regression function value estimator under Lipschitz class
#'
#' @inheritParams bias_Lip
#' @inheritParams var_Lip
#' @param TE logical specifying whether there are treatment and control groups.
#' @param d a vector of indicator variables specifying treatment and control group status;
#' relevant only when \code{TE = TRUE}.
#' @param alpha determines confidence level \code{1 - alpha}.
#' @param bw.eq if \code{TRUE}, the same bandwidths are used for estimators for treatment and control groups;
#' relevant only when \code{TE = TRUE}.
#' @param p Hölder exponent, either \code{p = 1} or \code{p = 2}. Default is \code{p = 1}.
#'
#' @return a list with the following components
#' \describe{
#'
#' \item{h.opt}{the optimal bandwidth; when \code{TE = TRUE}, return two bandwidths
#' for estimators for treatment and control groups.}
#'
#' \item{hl.opt}{the optimal half-length when the optimal bandwidth(s) is used.}
#'
#' \item{b.opt}{the bias corresponding to the optimal bandwidth \code{h.opt}.}
#'
#' \item{sd.opt}{the standard deviation corresponding to the optimal bandwidth
#' \code{h.opt}.}
#' }
#' @export
bw_opt <- function(y, x, t, TE = FALSE, d = NULL, M, kern, alpha, bw.eq = TRUE, p = 1){

  if(TE == TRUE){

    y.1 <- y[d == 1]
    y.0 <- y[d == 0]
    x.1 <- x[d == 1]
    x.0 <- x[d == 0]

    if(p == 1){
      b.fun <- function(h.1, h.0){
        bias_Lip(x.1, t, M, kern, h.1) + bias_Lip(x.0, t, M, kern, h.0)
      }
      sd.fun <- function(h.1, h.0){
        sqrt(var_Lip(y.1, x.1, t, kern, h.1) + var_Lip(y.0, x.0, t, kern, h.0))
      }
    }else if(p == 2){
      b.fun <- function(h.1, h.0){
        bias_Hol(x.1, t, M, kern, h.1) + bias_Hol(x.0, t, M, kern, h.0)
      }
      sd.fun <- function(h.1, h.0){
        sqrt(var_Hol(y.1, x.1, t, kern, h.1) + var_Hol(y.0, x.0, t, kern, h.0))
      }
    }

    obj <- function(h){

      h.1 <- abs(h[1]) # optim() might evaluate negative bandwidths
      h.0 <- abs(h[2])
      c <- stats::qnorm(1 - alpha) / 2
      return(b.fun(h.1, h.0) + c * sd.fun(h.1, h.0))
    }

    if(bw.eq == FALSE){

      h.1.init <- max(abs(x.1 - t)) / 2
      h.0.init <- max(abs(x.0 - t)) / 2

      opt.res <- stats::optim(c(h.1.init, h.0.init), obj, control = list(reltol = 1e-4))
      h.opt <- abs(opt.res$par) # optim() might evaluate negative bandwidths
      hl.opt <- opt.res$value  # half-length
    }else{

      obj.eq <- function(h){
        obj(c(h, h))
      }

      h.min <- max(sort(unique(abs(x.1 - t)))[2], sort(unique(abs(x.0 - t)))[2])
      h.max <- max(abs(x.1 - t), abs(x.0 - t))

      opt.res <- stats::optimize(obj.eq, c(h.min, h.max), tol = .Machine$double.eps^0.25)
      h.opt <- rep(opt.res$minimum, 2)
      hl.opt <- opt.res$objective
    }
    b.opt <- b.fun(h.opt[1], h.opt[2])
    sd.opt <- sd.fun(h.opt[1], h.opt[2])

  }else{

    try(if(p == 2) stop("use w_get_Hol()"))

    b.fun <- function(h){
      bias_Lip(x, t, M, kern, h)
    }
    sd.fun <- function(h){
      sqrt(var_Lip(y, x, t, kern, h))
    }

    obj.1 <- function(h){

      c <- stats::qnorm(1 - alpha) / 2
      return(b.fun(h) + c * sd.fun(h))
    }

    h.min <- sort(unique(abs(x - t)))[2]
    h.max <- abs(x - t)

    opt.res <- stats::optimize(obj.1, c(h.min, h.max), tol = .Machine$double.eps^0.25)
    h.opt <- opt.res$minimum
    hl.opt <- opt.res$objective
    b.opt <- b.fun(h.opt)
    sd.opt <- sd.fun(h.opt)
  }

  res <- list(h.opt = h.opt, hl.opt = hl.opt, b.opt = b.opt, sd.opt = sd.opt)
  return(res)
}


#' Optimal bandwidths under Lipschitz class with supplied conditional variance
#'
#' Calculates the optimal bandwidths for regression function value estimator under Lipschitz class
#' when the conditional variances are known.
#'
#' @inheritParams bias_Lip
#' @inheritParams var_Lip
#' @param TE logical specifying whether there are treatment and control groups.
#' @param d a vector of indicator variables specifying treatment and control group status;
#' relevant only when \code{TE = TRUE}.
#' @param c quantile value corresponding to \code{stats::qnorm(1 - alpha) / 2} given \code{alpha}
#' @param bw.eq if \code{TRUE}, the same bandwidths are used for estimators for treatment and control groups;
#' relevant only when \code{TE = TRUE}.
#' @param c.sd a vector of supplied conditional standard deviations
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
bw_Lip_supp <- function(c.sd, x, t, TE = FALSE, d = NULL, M, kern, c, bw.eq = TRUE){

  if(TE == TRUE){

    x.1 <- x[d == 1]
    x.0 <- x[d == 0]

    c.sd.1 <- c.sd[d == 1]
    c.sd.0 <- c.sd[d == 0]

    obj <- function(h){

      h.1 <- abs(h[1]) # optim() might evaluate negative bandwidths
      h.0 <- abs(h[2])
      bias <- bias_Lip(x.1, t, M, kern, h.1) + bias_Lip(x.0, t, M, kern, h.0)
      sd <- sqrt(var_Lip_resid(x.1, t, kern, h.1, c.sd.1) + var_Lip_resid(x.0, t, kern, h.0, c.sd.0))
      return(bias + c * sd)
    }

    if(bw.eq == FALSE){

      h.1.init <- max(abs(x.1 - t)) / 2
      h.0.init <- max(abs(x.0 - t)) / 2

      opt.res <- stats::optim(c(h.1.init, h.0.init), obj, control = list(reltol = 1e-4))
      h.opt <- abs(opt.res$par) # optim() might evaluate negative bandwidths
      hl.opt <- opt.res$value  # half-length
    }else{

      obj.eq <- function(h){
        obj(c(h, h))
      }

      h.min <- max(sort(unique(abs(x.1 - t)))[2], sort(unique(abs(x.0 - t)))[2])
      h.max <- max(abs(x.1 - t), abs(x.0 - t))

      opt.res <- stats::optimize(obj.eq, c(h.min, h.max), tol = .Machine$double.eps^0.25)
      h.opt <- rep(opt.res$minimum, 2)
      hl.opt <- opt.res$objective
    }

  }else{

    obj.1 <- function(h){

      bias <- bias_Lip(x, t, M, kern, h)
      sd <- sqrt(var_Lip_resid(x, t, kern, h, c.sd))
      return(bias + c * sd)
    }

    h.min <- sort(unique(abs(x - t)))[2]
    h.max <- abs(x - t)

    opt.res <- stats::optimize(obj.1, c(h.min, h.max), tol = .Machine$double.eps^0.25)
    h.opt <- opt.res$minimum
    hl.opt <- opt.res$objective
  }

  res <- list(h.opt = h.opt, hl.opt = hl.opt)
  return(res)
}
