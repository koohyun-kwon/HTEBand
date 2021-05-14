#' Optimal weights
#'
#' Solve for optimal weigths by the quantile matching method
#'
#' @param method character string used to define the parameter of interest and the function space.
#' Possible choices are
#' \describe{
#' \item{"reg.Hol"}{Regression function values under Hölder space}
#' \item{"reg.Lip"}{Regression function values under Lipschitz space}
#' \item{"TE.Lip"}{CATE under Lipschitz space}
#' \item{"TE.Lip.eqbw"}{CATE under Lipschitz space using the same bandwidths for the treatment
#' and the control groups}
#' }
#' @param C.vec smoothness parameter for the function space; possibly a vector.
#' @param y dependent variable; possibly a matrix with \code{nrow(y)} equals the number of observations.
#' @param x a vector of independent variable
#' @param d a vector of treatment indicator; it can be arbitrarily specified when not used.
#' @param eval a vector of indices
#' @param root.robust if \code{TRUE}, the fuction conducts diagnostic test whether optimizaiton worked well;
#' default is \code{root.robust = FALSE}.
#' @param ng the number of grids of quantile values over which the diagnostic test would be peformed
#' if \code{root.robust = TRUE}; default is \code{ng = 10}.
#' @inheritParams sup_quant_sim
#'
#' @return a list of the following components
#' \describe{
#' \item{w.1}{a 3-dim array of weigths corresponding to the treated observations;
#' the second and third dimension correponds to \code{ncol(y)} and \code{length(eval)}}
#' \item{w.0}{a 3-dim array of weigths corresponding to the control observations;
#' the second and third dimension correponds to \code{ncol(y)} and \code{length(eval)}.
#' If there is no control observations, returns array of 0.}
#' \item{c.root}{the quantile value corresponding to the optimal weigths.}
#' \item{increasing}{boolean value testing whether the objective function is increasing or not;
#' provided only if \code{root.robust = TRUE}}
#' \item{opt.grid}{a \code{data.frame} object containing evaluations of the objective function
#' over grid values of quantiles; provided only if \code{root.robust = TRUE}}
#' }
#' @export
opt_w <- function(method, C.vec, y, x, d, eval, T.grad.mat, level,
                  deg, kern, loo, M, seed = NULL, useloop = TRUE,
                  root.robust = FALSE, ng = 10){

  n.T <- length(eval)
  T.grad.mat <- v_to_m(T.grad.mat)
  k <- ncol(T.grad.mat)

  # Part 1: Setting initial variables

  if(method == "reg.Hol"){

    kern.reg <- "triangular"
    se.initial <- "EHW"
    se.method <- "nn"
    J <- 3
    C <- C.vec[1]
  }else if(method %in% c("reg.Lip", "TE.Lip", "TE.Lip.eqbw")){

    kern.reg <- "triangle"
    se.method <- "resid"
    C <- C.vec[1]
  }

  # Part 2: Calculating residuals

  if(method %in% c("reg.Hol", "reg.Lip")) d <- rep(1, length(y))

  y.1 <- v_to_m(y[d == 1])
  y.0 <-
    if(sum(d == 0) > 0){
      v_to_m(y[d == 0])
    }else{
      v_to_m(0)
    }

  x.1 <- x[d == 1]
  x.0 <-
    if(sum(d == 0) > 0){
      x[d == 0]
    }else{
      0
    }

  n.1 <- length(x.1)
  n.0 <- length(x.0)

  resid.1 <- matrix(0, nrow = n.1, ncol = k)
  resid.0 <- matrix(0, nrow = n.0, ncol = k)

  for(j in 1:k){

    resid.1[, j] <- eps_hat(y.1[, j], x.1, deg, kern, loo)
    if(sum(d == 0) > 0) resid.0[, j] <- eps_hat(y.0[, j], x.0, deg, kern, loo)
  }

  # Part 3: Defining objective function

  if(method %in% c("TE.Lip", "TE.Lip.eqbw")){

    # match orders of (y, x, d) and residuals
    y <- rbind(y.1, y.0)
    x <- c(x.1, x.0)
    d <- c(rep(1, n.1), rep(0, n.0))
    resid.1 <- as.vector(resid.1)
    resid.0 <- as.vector(resid.0)
    resid <- c(resid.1, resid.0)
  }

  eq <- function(c){

    level.int <- stats::pnorm(2 * c)

    if(method %in% c("reg.Hol", "reg.Lip", "TE.Lip", "TE.Lip.eqbw")){

      w.res <-
        if(method == "reg.Hol"){
          w_get_Hol(y, x, eval, C, level.int, kern.reg, se.initial, se.method, J)$w.mat
        }else if(method == "reg.Lip"){
          w_get_Lip(y, x, eval, C, level.int, kern = kern.reg, deg = deg, loo = loo,
                    se.method = se.method, resid = resid.1)
        }else if(method == "TE.Lip"){
          w_get_Lip(y, x, eval, C, level.int, TE = TRUE, d = d, kern = kern.reg,
                    bw.eq = FALSE, deg = deg, loo = loo, se.method = se.method, resid = resid)
        }else if(method == "TE.Lip.eqbw"){
          w_get_Lip(y, x, eval, C, level.int, TE = TRUE, d = d, kern = kern.reg,
                    bw.eq = TRUE, deg = deg, loo = loo, se.method = se.method, resid = resid)
        }

      w.1 <-
        if(sum(d == 0) > 0){
          w.res$w.mat.1
        }else{
          w.res
        }

      w.0 <-
        if(sum(d == 0) > 0){
          w.res$w.mat.0
        }else{
          0
        }

      w.1 <- array(w.1, dim = c(nrow(y.1), k, n.T))
      w.0 <- array(w.0, dim = c(nrow(y.0), k, n.T))

      q.sim <- sup_quant_sim(y.1, y.0, x.1, x.0, w.1, w.0, rep(1, n.T),
                             level, deg, kern, loo, M, seed, useloop, resid.1, resid.0)
    }

    eq.res <- list(val = c - q.sim, w.1 = w.1, w.0 = w.0)

    return(eq.res)
  }

  eq.val <- function(c){
    eq(c)$val
  }

  # Part 4: Optimization

  c.min <- stats::qnorm(level - 0.01)
  c.max <- stats::qnorm(1 - 0.1^3)

  root.res <- stats::uniroot(eq.val, interval = c(c.min, c.max), extendInt = "upX", tol = 1e-3)
  c.root <- root.res$root

  w.1 <- eq(c.root)$w.1
  w.0 <- eq(c.root)$w.0

  # Part 5: Optional step of robustness check for the optimization result

  if(root.robust == TRUE){

    c.grid <- seq(from = c.min, to = c.max, length.out = ng)
    obj.val <- numeric(ng)

    for(i in 1:ng) obj.val[i] = eq.val(c.grid[i])

    is.inc <- !is.unsorted(obj.val)
    opt.grid <- data.frame(c.grid = c.grid, obj.val = obj.val)

    res <- list(w.1 = w.1, w.0 = w.0, c.root = c.root, increasing = is.inc, opt.grid = opt.grid)

  }else{

    res <- list(w.1 = w.1, w.0 = w.0, c.root = c.root)
  }


  return(res)

}
