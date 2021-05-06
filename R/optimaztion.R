#' Optimal weights
#'
#' Solve for optimal weigths by the quantile matching method
#'
#' @param method character string used to define the parameter of interest and the function space.
#' Currently, only \code{"reg.Hol"} is supported.
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
  }else if(method == "reg.Lip" | method == "TE.Lip" | method == "TE.Lip.eqbw"){

    kern.reg <- "triangle"
    se.method <- "resid"
    C <- C.vec[1]
  }

  # Part 2: Calculating residuals

  if(method == "reg.Hol" | method == "reg.Lip"){

    y.1 <- v_to_m(y)
    y.0 <- v_to_m(0)
    n.1 <- length(y.1) / k
    n.0 <- length(y.0) / k

    resid.1 <- matrix(0, nrow = n.1, ncol = k)
    resid.0 <- matrix(0, nrow = n.0, ncol = k)

    for(j in 1:k){

      resid.1[, j] <- eps_hat(y.1[, j], x, deg, kern, loo)
      resid.0[, j] <- 0
    }
  }

  # Part 3: Defining objective function

  eq <- function(c){

    level.int <- stats::pnorm(2 * c)

    if(method == "reg.Hol" | method == "reg.Lip"){

      w.1 <-
        if(method == "reg.Hol"){
          w_get_Hol(y, x, eval, C, level.int, kern.reg, se.initial, se.method, J)$w.mat
        }else if(method == "reg.Lip"){
          w_get_Lip(y, x, eval, C, level.int, kern = kern.reg, deg = deg, loo = loo,
                    se.method = se.method)
        }

      w.1 <- array(w.1, dim = c(length(y), 1, n.T))
      w.0 <- array(0, dim = c(1, 1, n.T))

      q.sim <- sup_quant_sim(y, 0, x, 0, w.1, w.0, rep(1, n.T),
                             level, deg, kern, loo, M, seed, useloop, resid.1, resid.0)
    }

    eq.res <- list(val = c - q.sim, w.1 = w.1, w.0 = w.0)

    return(eq.res)
  }

  eq.val <- function(c){
    eq(c)$val
  }

  # Part 4: Optimization

  c.min <- stats::qnorm(level)
  c.max <- stats::qnorm(1 - (1 - level)/(2 * n.T)) + 1  # Add 1 to take into account some numerical errors

  root.res <- stats::uniroot(eq.val, interval = c(c.min, c.max))
  c.root <- root.res$root

  w.1 <- eq(c.root)$w.1
  w.0 <- eq(c.root)$w.0

  # Part 5: Optional step of robustness check for the optimization result

  if(root.robust){

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
