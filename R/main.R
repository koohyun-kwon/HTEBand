#' Confidence band for nonparametric regression
#'
#' Constructs a minimax optimal honest confidence band for nonparametric regression function values
#' under a specified function class.
#'
#' Currently, our method is guaranteed to work well for interior points. So setting \code{q.int > 0}
#' is recommended.
#'
#' @param y vector of dependent variables.
#' @param x vector of regressors.
#' @param C bound on the first or the second derivative, depending on \code{fclass}.
#' @param level coverage probability.
#' @param fclass function class specification, \code{"L"} for the Lipschitz class, and
#' \code{"H"} for the Hölder class.
#' @param n.eval number of grid points to use when constructing confidence band; default is
#' \code{n.eval = length(x) / 5}.
#' @param eval grid points to use when constructing confidence band; if it is not specified,
#' \code{eval} is automatically determined by \code{n.eval} and \code{q.int}.
#' @param q.int parameter determining the distance from the boundary of the support of \code{x}.
#' Confidence band is formed for points between the \code{q.int} quantile and the \code{1 - q.int} quantile
#' of \code{x}. \code{q.int} should be between 0 and 0.5. The default is \code{q.int = 0.025}.
#' @param n.sim number of bootstrap samples to use when calculating quantiles of suprema of empirical processes;
#' the default is \code{n.sim = 10^3}
#' @param kern type of kernel used in estimation; currently, only the triangle kernel (\code{kern = "triangle"})
#' is supported.
#' @param deg degree of local polynomial estimator used in the first-stage variance estimation;
#' the default is \code{deg = 1}.
#' @param loo logical indicating whether the leave-one-out procedure would be used in the first-stage variance estimation;
#' the default is \code{loo = FALSE}; see \code{\link[locpol]{locPolSmootherC}} for details.
#' @param seed seed number for bootstrap random sample generation. The default is \code{seed = NULL}.
#' @inheritParams cb_const
#'
#' @return  a data frame containing index set and corresponding confidence band values,
#' or if \code{root.robust = TRUE}, a list containing the data frame as well as
#' \code{increasing} and \code{opt.grid}; see \code{\link{opt_w}}.
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 500)
#' y <- x^2 + rnorm(500, 0, 1/4)
#' NpregBand(y, x, 2, 0.95, "L", n.eval = 25)
NpregBand <- function(y, x, C, level, fclass = c("L", "H"), n.eval = length(x) / 5, eval = NULL, q.int = 0.025,
                      n.sim = 10^3, kern = "triangle", deg = 1, loo = FALSE, seed = NULL,
                      root.robust = FALSE, ng = 10, x.out = NULL){

  fclass <- match.arg(fclass)
  method <-
    if(fclass == "L"){
      "reg.Lip"
    }else if(fclass == "H"){
      "reg.Hol"
    }

  if(is.null(eval)){
    eval.min <- stats::quantile(x, q.int)
    eval.max <- stats::quantile(x, 1 - q.int)
    eval <- seq(from = eval.min, to = eval.max, length.out = n.eval)
  }

  n.T <- length(eval)
  T.grad.mat <- rep(1, n.T)

  cb.res <- cb_const(method, C, y, x, 0, eval, T.grad.mat, level,
                     deg, kern, loo, n.sim, seed, useloop = TRUE,
                     root.robust, ng, x.out)

  return(cb.res)
}


#' Confidence band for conditional average treatment effects
#'
#' Constructs a minimax optimal honest confidence band for conditional average treatment effects
#' with respect to a regressor under a specified function class.
#'
#' Currently, our method is guaranteed to work well for interior points. So setting \code{q.int > 0}
#' is recommended.
#'
#' @inheritParams NpregBand
#' @param d vector of treatment indicators.
#' @param fclass function class specification; currently only \code{"L"} for the Lipschitz class is supported.
#' @param h.eq logical indicating whether the same bandwidths are used for the treated and the control groups;
#' the default is \code{h.eq = FALSE}.
#' @param n.eval number of grid points to use when constructing confidence band; default is
#' \code{min(sum(d == 1), sum(d == 0)) / 5}.
#'
#' @return a data frame containing index set and corresponding confidence band values,
#' or if \code{root.robust = TRUE}, a list containing the data frame as well as
#' \code{increasing} and \code{opt.grid}; see \code{\link{opt_w}}.
#' @export
#'
#' @examples
#' x <- rep(seq(-1, 1, length.out = 500), each = 2)
#' d <- rep(c(0, 1), 500)
#' y <- d * x^2 + (1 - d) * x + rnorm(500, 0, 1/4)
#' CATEBand(y, x, d, 2, 0.95, n.eval = 25)
CATEBand <- function(y, x, d, C, level, fclass = "L", h.eq = FALSE, n.eval = min(sum(d == 1), sum(d == 0)) / 5,
                     eval = NULL, q.int = 0.025, n.sim = 10^3, kern = "triangle", deg = 1, loo = FALSE, seed = NULL,
                     root.robust = FALSE, ng = 10, x.out = NULL){

  method <-
    if(h.eq == TRUE){
      "TE.Lip.eqbw"
    }else{
      "TE.Lip"
    }

  if(is.null(eval)){
    x.1 <- x[d == 1]
    x.0 <- x[d == 0]
    eval.min <- max(stats::quantile(x.1, q.int), stats::quantile(x.0, q.int))
    eval.max <- min(stats::quantile(x.1, 1 - q.int), stats::quantile(x.0, 1 - q.int))
    eval <- seq(from = eval.min, to = eval.max, length.out = n.eval)
  }

  n.T <- length(eval)
  T.grad.mat <- rep(1, n.T)

  cb.res <- cb_const(method, C, y, x, d, eval, T.grad.mat, level,
                     deg, kern, loo, n.sim, seed, useloop = TRUE,
                     root.robust, ng, x.out)

  return(cb.res)
}
