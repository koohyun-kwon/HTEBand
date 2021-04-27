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
avar <- function(w.1, w.0, omega.1, omega.0, T.grad){

  k <- length(T.grad)

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

  dmat.1 <- matrix(0, k, k)
  dmat.0 <- dmat.1

  for(j1 in 1:k){
    for(j2 in j1:k){

      dmat.1[j1, j2] <- T.grad[j1] * T.grad[j2] * sum(w.1[, j1] * w.1[, j2] * omega.1[, j1, j2])
      dmat.0[j1, j2] <- T.grad[j1] * T.grad[j2] * sum(w.0[, j1] * w.0[, j2] * omega.0[, j1, j2])

      if(j1 != j2){

        dmat.1[j2, j1] <- dmat.1[j1, j2]
        dmat.0[j2, j1] <- dmat.0[j1, j2]
      }
    }

  }

  avar <- sum(dmat.1) + sum(dmat.0)
  dnmnt <- sqrt(avar)

  return(dnmnt)
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

# stud_err_sim <- function(w.1, w.0, y.1, y.0, x.1, x.0, T.grad, deg, kern, loo){
#
#
# }

