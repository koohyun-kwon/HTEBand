#' Studentized error variable
#'
#' Used to calculate the quantile of the supremum of these error variables over different values of indices.
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

    if(!is.matrix(w.1)){

      w.1 <- matrix(w.1, length(w.1), 1)
    }

    if(!is.matrix(w.0)){

      w.0 <- matrix(w.0, length(w.0), 1)
    }

    if(!is.matrix(resid.1)){

      resid.1 <- matrix(resid.1, length(resid.1), 1)
    }

    if(!is.matrix(resid.0)){

      resid.0 <- matrix(resid.0, length(resid.0), 1)
    }

    if(length(dim(omega.1)) != 3){

      if(length(dim(omega.1)) == 0){

        omega.1 <- array(omega.1, dim = c(length(omega.1), k, k))
      }else{

        stop("omega.1 should be either a 3-dim array or a vector.")
      }
    }

    if(length(dim(omega.0)) != 3){

      if(length(dim(omega.0)) == 0){

        omega.0 <- array(omega.0, dim = c(length(omega.0), k, k))
      }else{

        stop("omega.0 should be either a 3-dim array or a vector.")
      }
    }
  }


  # Numerator

  nmrt <- sum(T.grad * diag(t(w.1) %*% resid.1 - t(w.0) %*% resid.0))

  # Denominator

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

  # Final result
  res <- nmrt / dnmnt

  return(res)
}
