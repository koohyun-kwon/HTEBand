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
#' \item{"TE.Lip"}{CATE under Hölder space}
#' \item{"TE.Lip.eqbw"}{CATE under Hölder space using the same bandwidths for the treatment
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
#' @param resid.1 residuals corresponding to the treated observations; it can be a vector
#' or a matrix if \code{ncol(y) > 1}. For treatment effect models, its \eqn{i}th component
#' corresponds to the \eqn{i}th component of \code{y[d == 1, ]}. If there are no
#' treatment and control groups, \code{resid.1} corresponds to residuals for the
#' entire observations. It can be specified to be \code{NULL}.
#' @param resid.0 residuals corresponding to the control observations; it can be a vector
#' or a matrix if \code{ncol(y) > 1}. Its \eqn{i}th component
#' corresponds to the \eqn{i}th component of \code{y[d == 0, ]}. It can be specified to be
#' \code{NULL}.
#' @param c.method method to calculate the optimal value of \eqn{c_n}, an element
#' in \code{c("root", "supp")}; default is \code{c.method = "root"}.
#' @param c.supp value of \eqn{c_n} to be used when \code{c.method = "supp"};
#' default is \code{c.supp = NULL}.
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
opt_w <- function(method, C.vec, y, x, d = NULL, eval, T.grad.mat, level,
                  deg, kern, M, seed = NULL, useloop = TRUE,
                  root.robust = FALSE, ng = 10, resid.1 = NULL, resid.0 = NULL,
                  var.reg = "npr", c.method = "root", c.supp = NULL){

  n.T <- length(eval)
  T.grad.mat <- v_to_m(T.grad.mat)
  k <- ncol(T.grad.mat)
  is.TE <- 0 %in% d & sum(d == 0) > 1
  y <- v_to_m(y)
  n <- nrow(y)
  if(!is.TE) d <- rep(1, n)

  y.1 <- v_to_m(y[d == 1, ])
  x.1 <- x[d == 1]

  if(is.TE){

    # Reordering (y, x, d): treated obs first, control obs later
    y.0 <- v_to_m(y[d == 0, ])
    y <- rbind(y.1, y.0)
    x.0 <- x[d == 0]
    x <- c(x.1, x.0)
    d <- c(rep(1, length(x.1)), rep(0, length(x.0)))
  }else{

    # Even if there is no control group, these values will be used as an input
    # in the quantile calculation function
    y.0 <- matrix(0, nrow = 1, ncol = k)
    x.0 <- 0
  }


  # Part 1: Setting initial variables
  if(method %in% c("reg.Hol", "TE.Hol", "TE.Hol.eqbw")){

    kern.reg <- "triangular"
    se.initial <- "EHW"
    se.method <- "nn"
    J <- 3
    C <- C.vec[1]
  }else if(method %in% c("reg.Lip", "TE.Lip", "TE.Lip.eqbw")){

    kern.reg <- "tri"
    se.method <- "resid"
    C <- C.vec[1]
  }

  # Part 2: Calculating residuals

  if(is.null(resid.1)){

    resid.all <- resid_calc(y, x, d, deg, kern, var.reg)
    resid.1 <- resid.all$resid.1
    resid.0 <- resid.all$resid.0

  }else{

    resid.1 <- v_to_m(resid.1)
    if(!is.null(resid.0)){
      resid.0 <- v_to_m(resid.0)
    }else{
      resid.0 <- matrix(0, nrow = 1, ncol = k)
    }
  }

  # Part 3: optimal weight function and simulated critical value

  eq <- function(c){

    level.int <- stats::pnorm(2 * c)

    if(method %in% c("reg.Hol", "reg.Lip", "TE.Lip", "TE.Lip.eqbw",
                     "TE.Hol", "TE.Hol.eqbw")){

      w.res <-
        if(method == "reg.Hol"){
          w_get_Hol(y, x, eval, C, level.int, kern.reg, se.initial, se.method, J)
        }else if(method == "reg.Lip"){
          w_get_Lip(y, x, eval, C, level.int, kern = kern.reg)
        }else if(method == "TE.Lip"){
          w_get_Lip(y, x, eval, C, level.int, TE = TRUE, d = d, kern = kern.reg,
                    bw.eq = FALSE)
        }else if(method == "TE.Lip.eqbw"){
          w_get_Lip(y, x, eval, C, level.int, TE = TRUE, d = d, kern = kern.reg,
                    bw.eq = TRUE)
        }else if(method == "TE.Hol"){
          w_get_Hol(y, x, eval, C, level, kern.reg, se.initial,
                    se.method, J, TE = TRUE, d = d, bw.eq = FALSE)
        }else if(method == "TE.Hol.eqbw"){
          w_get_Hol(y, x, eval, C, level, kern.reg, se.initial,
                    se.method, J, TE = TRUE, d = d, bw.eq = TRUE)
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
                             level, deg, kern, M, seed, useloop, resid.1, resid.0)
    }

    # eq.res <- list(val = c - q.sim, w.1 = w.1, w.0 = w.0)
    eq.res <- list(q.sim = q.sim, w.1 = w.1, w.0 = w.0)

    return(eq.res)
  }

  # Part 4: Weight function calculation

  if(c.method == "root"){

    eq.val <- function(c){
      c - eq(c)$q.sim
    }

    c.min <- stats::qnorm(level - 0.01)
    c.max <- stats::qnorm(1 - 0.1^3)

    root.res <- stats::uniroot(eq.val, interval = c(c.min, c.max), extendInt = "upX", tol = 1e-3)
    c.root <- root.res$root

  }else if(c.method %in% c("supp", "supp2")){

    try(if(is.null(c.supp)) stop("Supply value of c.supp"))
    c.root <- c.supp
  }

  eq.res <- eq(c.root)

  w.1 <- eq.res$w.1
  w.0 <- eq.res$w.0
  q.sim <- eq.res$q.sim

  # Part 5: Optional step of robustness check for the optimization result

  if(root.robust == TRUE & c.method == "root"){

    c.grid <- seq(from = c.min, to = c.max, length.out = ng)
    obj.val <- numeric(ng)

    for(i in 1:ng) obj.val[i] = eq.val(c.grid[i])

    is.inc <- !is.unsorted(obj.val)
    opt.grid <- data.frame(c.grid = c.grid, obj.val = obj.val)

    res <- list(w.1 = w.1, w.0 = w.0, c.root = c.root, q.sim = q.sim,
                increasing = is.inc, opt.grid = opt.grid)

  }else{

    res <- list(w.1 = w.1, w.0 = w.0, c.root = c.root, q.sim = q.sim)
  }


  return(res)

}
