#' Optimal confidence band
#'
#' Calculate confidence band by searching for the optimal coverage probability
#' for individual confidence intervals.
#'
#' The confidence band is calculated over the points in \code{eval} and interpolated over
#' \code{x.out} using linear approximation.
#'
#' @inheritParams opt_w
#' @param x.out the grid of points the confidence band is evaluated; see details.
#'
#' @return a data frame containing index set, corresponding confidence band values, and
#' bandwidths used for the treatment and the control groups (for models without treatment and control groups,
#' the same value of bandwidths will be returned). If \code{root.robust = TRUE},
#' the function returns a list containing the data frame without bandwidths data
#'  as well as \code{increasing} and \code{opt.grid}; see \code{\link{opt_w}}.
#' @export
cb_const <- function(method, C.vec, y, x, d, eval, T.grad.mat, level,
                     deg, kern, M, var.reg = "npr", seed = NULL, useloop = TRUE,
                     root.robust = FALSE, ng = 10, x.out = NULL,
                     c.method = "root"){

  n.T <- length(eval)
  cb.grid <- matrix(NA, nrow = n.T, ncol = 4)
  n <- length(x)

  if(method %in% c("reg.Hol", "TE.Hol", "TE.Hol.eqbw")){

    kern.reg <- "triangular"
    se.initial <- "EHW"
    se.method <- "nn"
    J <- 3
    C <- C.vec[1]
    p <- 2
  }else if(method %in% c("reg.Lip", "TE.Lip", "TE.Lip.eqbw")){

    kern.reg <- "tri"
    se.method <- "resid"
    C <- C.vec[1]
    p <- 1
  }

  resid.res <- resid_calc(y, x, d, deg, kern, var.reg)
  resid.1 <- resid.res$resid.1
  resid.0 <- resid.res$resid.0

  if(c.method == "root"){
    c.supp = NULL
  }else if(c.method == "supp"){
    c.supp = opt_cn(n, p)
  }

  opt.res <- opt_w(method, C.vec, y, x, d, eval, T.grad.mat, level,
                   deg, kern, M, seed, useloop,
                   root.robust, ng, resid.1, resid.0, c.method = c.method,
                   c.supp = c.supp)

  c.opt <- opt.res$c.root
  ci.level <- stats::pnorm(2 * c.opt)

  if(c.method == "root"){
    ci.cv = NULL
  }else if(c.method == "supp"){
    ci.cv <- opt.res$q.sim
  }

  if(root.robust & c.method == "root"){
    increasing <- opt.res$increasing
    opt.grid <- opt.res$opt.grid
  }

  for(t in 1:n.T){

    cb.grid[t, ] <-
      if(method == "reg.Hol"){
        ci_reg_Hol(y, x, eval[t], C, ci.level, kern.reg, se.initial, se.method, J,
                   cv = ci.cv)
      }else if(method == "reg.Lip"){
        ci_reg_Lip(y, x, eval[t], C, ci.level, kern = kern.reg,
                   se.method = se.method, cv = ci.cv)
      }else if(method == "TE.Lip"){
        ci_reg_Lip(y, x, eval[t], C, ci.level, TE = TRUE, d = d, kern = kern.reg,
                   bw.eq = FALSE, se.method = se.method, cv = ci.cv)
      }else if(method == "TE.Lip.eqbw"){
        ci_reg_Lip(y, x, eval[t], C, ci.level, TE = TRUE, d = d, kern = kern.reg,
                   bw.eq = TRUE, se.method = se.method, cv = ci.cv)
      }else if(method == "TE.Hol"){
        ci_reg_Hol(y, x, eval[t], C, ci.level, TE = TRUE, d = d, kern = kern.reg,
                   bw.eq = FALSE, se.method = se.method, cv = ci.cv)
      }else if(method == "TE.Hol.eqbw"){
        ci_reg_Hol(y, x, eval[t], C, ci.level, TE = TRUE, d = d, kern = kern.reg,
                   bw.eq = TRUE, se.method = se.method, cv = ci.cv)
      }
  }

  if(is.null(x.out)){
    cb.data <- data.frame(eval = eval, cb.lower = cb.grid[, 1], cb.upper = cb.grid[, 2],
                          h.t = cb.grid[, 3], h.c = cb.grid[, 4])
  }else{
    cb.l <- stats::approx(eval, cb.grid[, 1], x.out)$y
    cb.u <- stats::approx(eval, cb.grid[, 2], x.out)$y
    cb.data <- data.frame(x.out = x.out, cb.lower = cb.l, cb.upper = cb.u)
  }

  res <-
    if(root.robust){
      list(cb.data = cb.data, increasing = increasing, opt.grid = opt.grid)
    }else{
      cb.data
    }

  return(res)
}
