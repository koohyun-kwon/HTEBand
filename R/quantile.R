#' Asymptotic variance
#'
#' Calculates the asymptotic variance expression given weights and conditional variance matrix.
#'
#' @param w.1 A \code{n_1} by \code{k} dimensional matrix of weight values corresponding to treated observations
#' @param w.0 A \code{n_0} by \code{k} dimensional matrix of weight values corresponding to control observations
#' @param omega.1 A \code{n_1} by \code{k} by \code{k} array of conditional variance corresponding to treated observations
#' @param omega.0 A \code{n_0} by \code{k} by \code{k} array of conditional variance corresponding to control observations
#' @param T.grad A \code{k} dimensional gradient vector of \code{T_t f}
#'
#' @return A scalar asymptotic variance value
#' @export
#'
#' @examples
#' avar(rep(1/50, 50), rep(1/50, 50), rep(1, 50), rep(1, 50), 1)
#' w.1 <- w.0 <- matrix(rep(1/50, 100), ncol = 2)
#' omega.1 <- omega.0 <- array(rep(c(1,0,0,1), each = 50), dim = c(50, 2, 2))
#' T.grad <- c(1, 1)
#' avar(w.1, w.0, omega.1, omega.0, T.grad)
avar <- function(w.1, w.0, omega.1, omega.0, T.grad){

  k <- length(T.grad)
  n.1 <- length(w.1) / k
  n.0 <- length(w.0) / k

  # handling vector inputs when k = 1

  if(k == 1){

    w.1 <- v_to_m(w.1)
    w.0 <- v_to_m(w.0)

    omega.1 <- v_to_3d(omega.1, "omega.1")$res.arr
    omega.0 <- v_to_3d(omega.0, "omega.0")$res.arr

    if(v_to_3d(omega.1, "omega.1")$err.stat == 1){
      stop(v_to_3d(omega.1, "omega.1")$err.msg)
    }

    if(v_to_3d(omega.0, "omega.0")$err.stat == 1){
      stop(v_to_3d(omega.0, "omega.0")$err.msg)
    }
  }

  T.grad.rep.1 <- matrix(rep(T.grad, each = n.1), n.1, k)
  T.grad.rep.0 <- matrix(rep(T.grad, each = n.0), n.0, k)
  T.grad.arr.1 <- mat_sq(T.grad.rep.1)
  T.grad.arr.0 <- mat_sq(T.grad.rep.0)

  w.arr.1 <- mat_sq(w.1)
  w.arr.0 <- mat_sq(w.0)

  res <- sum(T.grad.arr.1 * w.arr.1 * omega.1) + sum(T.grad.arr.0 * w.arr.0 * omega.0)

  return(sqrt(res))
}

#' Studentized error variable
#'
#' Calculates the studentized error term given weights, residuals, and conditional variance matrix.
#'
#' @param w.1 A \code{n_1} by \code{k} dimensional matrix of weight values corresponding to treated observations
#' @param w.0 A \code{n_0} by \code{k} dimensional matrix of weight values corresponding to control observations
#' @param resid.1 A \code{n_1} by \code{k} dimensional matrix of residuals corresponding to treated observations
#' @param resid.0 A \code{n_0} by \code{k} dimensional matrix of residuals corresponding to control observations
#' @param omega.1 A \code{n_1} by \code{k} by \code{k} array of conditional variance corresponding to treated observations
#' @param omega.0 A \code{n_0} by \code{k} by \code{k} array of conditional variance corresponding to control observations
#' @param T.grad A \code{k} dimensional gradient vector of \code{T_t f}
#'
#' @return A studentized scalar estimation error term
#' @export
#'
#' @examples
#' stud_err(rep(1/50, 50), rep(1/50, 50), stats::rnorm(50),
#' stats::rnorm(50), rep(1, 50), rep(1, 50), 1)
#' w.1 <- w.0 <- matrix(rep(1/50, 100), ncol = 2)
#' resid.1 <- resid.0 <- matrix(stats::rnorm(100), ncol = 2)
#' omega.1 <- omega.0 <- array(rep(c(1,0,0,1), each = 50), dim = c(50, 2, 2))
#' T.grad <- c(1, 1)
#' stud_err(w.1, w.0, resid.1, resid.0, omega.1, omega.0, T.grad)
stud_err <- function(w.1, w.0, resid.1, resid.0, omega.1, omega.0, T.grad){

  k <- length(T.grad)

  # handling vector inputs when k = 1

  if(k == 1){

    w.1 <- v_to_m(w.1)
    w.0 <- v_to_m(w.0)
    resid.1 <- v_to_m(resid.1)
    resid.0 <- v_to_m(resid.0)

    omega.1 <- v_to_3d(omega.1, "omega.1")$res.arr
    omega.0 <- v_to_3d(omega.0, "omega.0")$res.arr

    if(v_to_3d(omega.1, "omega.1")$err.stat == 1){
      stop(v_to_3d(omega.1, "omega.1")$err.msg)
    }

    if(v_to_3d(omega.0, "omega.0")$err.stat == 1){
      stop(v_to_3d(omega.0, "omega.0")$err.msg)
    }
  }


  # Numerator

  nmrt <- sum(T.grad * diag(t(w.1) %*% resid.1 - t(w.0) %*% resid.0))

  # Denominator

  dnmnt <- avar(w.1, w.0, omega.1, omega.0, T.grad)

  # Final result
  res <- nmrt / dnmnt

  return(res)
}

# stud_err_sim <- function(y.1, y.0, x.1, x.0, w.1, w.0, T.grad, deg, kern, loo, M, seed = 1){
#
#   k <- length(T.grad)
#   n.1 <- length(y.1) / k
#   n.0 <- length(y.0) / k
#
#   set.seed(seed)
#   z.1 <- matrix(stats::rnorm(n.1 * K * M), n.1, k * M)
#   z.0 <- matrix(stats::rnorm(n.0 * K * M), n.0, k * M)
#
#   w.1.rep <- matrix(rep(w.1, M), n.1, k * M)
#   w.0.rep <- matrix(rep(w.0, M), n.0, k * M)
#   T.grad.1.rep <- matrix(rep(matrix(rep(T.grad, each = n.1), ncol = k), M), n.1, k * M)
#   T.grad.0.rep <- matrix(rep(matrix(rep(T.grad, each = n.0), ncol = k), M), n.0, k * M)
#
#   resid.1 <- matrix(0, nrow = n.1, ncol = k)
#   resid.0 <- matrix(0, nrow = n.0, ncol = k)
#   y.1 <- v_to_m(y.1)
#   y.0 <- v_to_m(y.0)
#
#   for(j in 1:k){
#
#     resid.1[, j] <- eps_hat(y.1[, j], x, deg, kern, loo)
#     resid.0[, j] <- eps_hat(y.0[, j], x, deg, kern, loo)
#   }
#
#   resid.1.rep <- matrix(rep(resid.1, M), n.1, k * M)
#   resid.0.rep <- matrix(rep(resid.0, M), n.0, k * M)
#
#   nmrt.pre.1 <- T.grad.1.rep * w.1.rep * resid.1.rep
#   nmrt.pre.0 <- T.grad.0.rep * w.0.rep * resid.0.rep
#
#   nmrt.1 <- colSums(matrix(colSums(nmrt.pre.1), k, M))
#   nmrt.0 <- colSums(matrix(colSums(nmrt.pre.0), k, M))
#   nmrt <- nmrt.1 - nmrt.0
#
#   omega.hat.1 <- mat
#
#
# }
#
