#' Optimal confidence band
#'
#' Calculate confidence band by searching for the optimal coverage probability
#' for individual confidence intervals.
#'
#' @inheritParams opt_w
#'
#' @return a data frame containing index set and corresponding confidence band values,
#' or if \code{root.robust = TRUE}, a list containing the data frame as well as
#' \code{increasing} and \code{opt.grid}; see \code{\link{opt_w}}.
#' @export
cb_const <- function(method, C.vec, y, x, d, eval, T.grad.mat, level,
                     deg, kern, loo, M, seed = NULL, useloop = TRUE,
                     root.robust = FALSE, ng = 10){

  n.T <- length(eval)
  cb.data <- data.frame(eval = eval, cb.lower = numeric(n.T), cb.upper = numeric(n.T))

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

  opt.res <- opt_w(method, C.vec, y, x, d, eval, T.grad.mat, level,
                   deg, kern, loo, M, seed, useloop,
                   root.robust, ng)
  c.opt <- opt.res$c.root

  if(root.robust){
    increasing <- opt.res$increasing
    opt.grid <- opt.res$opt.grid
  }

  ci.level <- stats::pnorm(2 * c.opt)

  for(t in 1:n.T){

    cb.data[t, 2:3] <-
      if(method == "reg.Hol"){
        ci_reg_Hol(y, x, eval[t], C, ci.level, kern.reg, se.initial, se.method, J)
      }else if(method == "reg.Lip"){
        ci_reg_Lip(y, x, eval[t], C, ci.level, kern = kern.reg,
                   deg = deg, loo = loo, se.method = se.method)
      }else if(method == "TE.Lip"){
        ci_reg_Lip(y, x, eval[t], C, ci.level, TE = TRUE, d = d, kern = kern.reg,
                   bw.eq = FALSE, deg = deg, loo = loo, se.method = se.method)
      }else if(method == "TE.Lip.bweq"){
        ci_reg_Lip(y, x, eval[t], C, ci.level, TE = TRUE, d = d, kern = kern.reg,
                   bw.eq = TRUE, deg = deg, loo = loo, se.method = se.method)
      }
  }

  res <-
    if(root.robust){
      list(cb.data = cb.data, increasing = increasing, opt.grid = opt.grid)
    }else{
      cb.data
    }

  return(res)
}
