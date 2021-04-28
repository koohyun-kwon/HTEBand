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

  return(res)
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

  dnmnt <- sqrt(avar(w.1, w.0, omega.1, omega.0, T.grad))

  # Final result
  res <- nmrt / dnmnt

  return(res)
}

#' Simulation draws of studentized error
#'
#' Produce a vector of simulation draws of studentized errors
#'
#' @param y.1 dependent variable for treated observation; possibly a matrix with \code{nrow} equals the sample size
#' @param y.0 dependent variable for control observation; possibly a matrix with \code{nrow} equals the sample size
#' @param x.1 a vector of regressor for treated observation
#' @param x.0 a vector of regressor for control observation
#' @inheritParams stud_err
#' @inheritParams eps_hat
#' @param z.1 a vector of simulated standard normal random variables of length \code{length(y.1) * M} for some \code{M}
#' @param z.0 a vector of simulated standard normal random variables of length \code{length(y.0) * M} for some \code{M}
#'
#' @return a list of following components
#' \describe{
#' \item{err.sim}{\code{M} dimensional vector of simulated studentized errors}
#' \item{nmrt}{numerator for err.sim}
#' \item{dnmnt}{denominator for err.sim}
#' }
#'
#' @export
#'
#' @examples
#' x.1 <- x.0 <- seq(from = -1, to = 1, length.out = 500)
#' y.1 <- x.1^2 + stats::rnorm(500, 0, 0.1)
#' y.0 <- x.0^2 + stats::rnorm(500, 0, 0.1)
#' w.1 <- w.0 <- rep(1/500, 500)
#' z.1 <- rnorm(500 * 500)
#' z.0 <- rnorm(500 * 500)
#' stud_err_sim(y.1, y.0, x.1, x.0, w.1, w.0, 1, 1, "triangle", TRUE, z.1, z.0)
stud_err_sim <- function(y.1, y.0, x.1, x.0, w.1, w.0, T.grad, deg, kern, loo, z.1, z.0){

  k <- length(T.grad)
  n.1 <- length(y.1) / k
  n.0 <- length(y.0) / k
  M <- length(z.1) / (k * n.1)

  z.1 <- matrix(z.1, n.1, k * M)
  z.0 <- matrix(z.0, n.0, k * M)

  w.1.rep <- matrix(rep(w.1, M), n.1, k * M)
  w.0.rep <- matrix(rep(w.0, M), n.0, k * M)
  T.grad.1.rep <- matrix(rep(matrix(rep(T.grad, each = n.1), ncol = k), M), n.1, k * M)
  T.grad.0.rep <- matrix(rep(matrix(rep(T.grad, each = n.0), ncol = k), M), n.0, k * M)

  resid.1 <- matrix(0, nrow = n.1, ncol = k)
  resid.0 <- matrix(0, nrow = n.0, ncol = k)
  y.1 <- v_to_m(y.1)
  y.0 <- v_to_m(y.0)

  for(j in 1:k){

    resid.1[, j] <- eps_hat(y.1[, j], x.1, deg, kern, loo)

    if(n.0 == 1){
      resid.0[, j] <- 0
    }else{
      resid.0[, j] <- eps_hat(y.0[, j], x.0, deg, kern, loo)
    }

  }

  resid.1.rep <- matrix(rep(resid.1, M), n.1, k * M)
  resid.0.rep <- matrix(rep(resid.0, M), n.0, k * M)

  nmrt.pre.1 <- T.grad.1.rep * w.1.rep * resid.1.rep * z.1
  nmrt.pre.0 <- T.grad.0.rep * w.0.rep * resid.0.rep * z.0

  nmrt.1 <- colSums(matrix(colSums(nmrt.pre.1), k, M))
  nmrt.0 <- colSums(matrix(colSums(nmrt.pre.0), k, M))
  nmrt <- nmrt.1 - nmrt.0

  omega.hat.1 <- mat_sq(resid.1)
  omega.hat.0 <- mat_sq(resid.0)

  dnmnt <- sqrt(avar(w.1, w.0, omega.hat.1, omega.hat.0, T.grad))

  return(list(err.sim = nmrt / dnmnt, nmrt = nmrt, dnmnt = dnmnt))
}

#' Quantile of supremum of a studentized process
#'
#' Calculates the quantile of supremum of the absolute value of a studentized process
#' indexed by a index set.
#'
#' @inheritParams stud_err_sim
#' @param w.1.arr A \code{n_1} by \code{k} by \code{n.T} dimensional array of weight values corresponding to treated observations
#' @param w.0.arr A \code{n_0} by \code{k} by \code{n.T} dimensional array of weight values corresponding to treated observations
#' @param T.grad.mat A \code{n.T} by \code{k} dimensional gradient matrix of \code{T_t f} for \code{t} = 1,..., \code{n.T}.
#' @param level level of quantile
#' @param useloop If \code{TRUE}, the function is implemented by \code{for} loop over \code{t} = 1,..., \code{n.T}.
#'
#' @return a scalar quantile value
#' @export
#'
#' @examples
#' n <- 500
#' x <- stats::runif(n, min = -1, max = 1)
#' y <- x + rnorm(n, 0, 1/4)
#' n.T <- 10
#' eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
#' w <- array(w_get_Hol(y, x, eval, 1, 0.95)$w.mat, dim = c(n, 1, n.T))
#' z <- rnorm(n * 100)
#' sup_quant_sim(y, 0, x, 0, w, array(rep(0, n.T), dim = c(1, 1, n.T)),
#' rep(1, n.T), 0.95, 1, "triangle", FALSE, z, rnorm(100), useloop = TRUE)
sup_quant_sim <- function(y.1, y.0, x.1, x.0, w.1.arr, w.0.arr, T.grad.mat, level,
                        deg, kern, loo, z.1, z.0, useloop = TRUE){

  T.grad.mat <- v_to_m(T.grad.mat)
  n.T <- nrow(T.grad.mat)
  k <- ncol(T.grad.mat)
  n.1 <- length(y.1) / k
  n.0 <- length(y.0) / k
  M <- length(z.1) / (k * n.1)

  if(useloop){

    max.val <- rep(0, M)
    for(t in 1:n.T){

      T.grad <- T.grad.mat[t, ]
      w.1 <- w.1.arr[, , t]
      w.0 <- w.0.arr[, , t]
      val.new <- abs(stud_err_sim(y.1, y.0, x.1, x.0, w.1, w.0, T.grad, deg, kern, loo, z.1, z.0)$err.sim)
      max.val <- pmax(max.val, val.new)
    }
  }

  return(stats::quantile(max.val, level))
}

