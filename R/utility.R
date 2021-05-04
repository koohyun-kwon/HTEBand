#' Vector to matrix
#'
#' Test whether an input is a vector, and if so, transform it into a matrix with \code{ncol = 1}.
#'
#' @param v a matrix or a vector
#'
#' @return a corresponding matrix
#' @export
#'
#' @examples
#' v_to_m(1:10)
#' v_to_m(matrix(1:10, ncol = 1))
v_to_m <- function(v){

  if(!is.matrix(v)){
     res <- matrix(v, length(v), 1)
  }else{
    res <- v
  }

  return(res)
}

#' Vector to 3-dim array
#'
#'  Test whether an input is a vector, and if so, transform it into a array with its second and third dimension lengths equal to 1.
#'  If an input is neither a vector nor a 3-dimensional array, returns \code{NULL} with an error message.
#'
#' @param v a vector or an array
#' @param obj.name a string used to display the error message
#'
#' @return a list containing the following components
#' \describe{
#'
#'   \item{\code{res.arr}}{Corresponding 3-dimensional array. If the input is not right, this contains \code{NULL}}
#'
#'   \item{\code{err.stat}}{Error status variable, taking 1 if the input it not right. Otherwise, this equals 0.}
#'
#'   \item{\code{err.msg}}{Error message to display if \code{err.stat = 1}. This contains an empty string if \code{err.stat = 0}.}
#'
#' }
#' @export
#'
#' @examples
#'
#' v_to_3d(1:10, "omega")
#' v_to_3d(array(1:10, dim = c(10, 1, 1)), "omega")
#' v_to_3d(matrix(1:10, ncol = 1), "omega")
v_to_3d <- function(v, obj.name){

  res.arr <- v
  err.stat <- 0
  err.msg <- ""

  if(length(dim(res.arr)) != 3){

    if(length(dim(res.arr)) == 0){

      res.arr <- array(res.arr, dim = c(length(res.arr), 1, 1))
    }else{

      res.arr <- NULL
      err.stat <- 1
      err.msg <- paste(obj.name, "should be a 3-dim array or a vector")
    }
  }

  res <- list(res.arr = res.arr, err.stat = err.stat, err.msg = err.msg)

  return(res)
}

#' Columnwise matrix self-multiplication
#'
#' Given a \code{n} by \code{k} matrix \code{mat}, computes a  \code{n} by \code{k} by \code{k} array
#' whose \code{(i, j1, j2)}th component is given by \code{mat[i, j1] * mat[i, j2]}.
#'
#' @param mat A matrix or a vector
#' @param arr.ret If \code{FALSE}, returns \code{n} by \code{k^2} matrix instead. Default is \code{TRUE}.
#'
#' @return A \code{n} by \code{k} by \code{k} array or a \code{n} by \code{k^2} matrix.
#' @export
#'
#' @examples
#' mat <- matrix(1:12, ncol = 3)
#' mat_sq(mat, TRUE)
#' mat_sq(1:10, FALSE)
mat_sq <- function(mat, arr.ret = TRUE){

  mat <- v_to_m(mat)
  k <- ncol(mat)
  n <- nrow(mat)

  mat.1 <- t(matrix(rep(t(mat), each = k), nrow = k^2))
  mat.2 <- matrix(rep(mat, k), ncol = k^2)

  if(arr.ret){
    res <- array(mat.1 * mat.2, dim = c(n, k, k))
  }else{
    res <- mat.1 * mat.2
  }

  return(res)
}

#' Kernel function
#'
#' Calculates the values of kernel function given a vector of regressors.
#'
#' @param x a vector of regressors.
#' @param t a scalar evaluation point.
#' @param h a positive bandwidth.
#' @param kern a string for kernel name; currently \code{"triangle"} is supported.
#'
#' @return a vector of kernel values with the same dimension as \code{x}.
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 10)
#' K_fun(x, 0, 0.5, "triangle")
K_fun <- function(x, t, h, kern = c("triangle")){

  if(is.na(h) | h <=0) stop("Invalid bandwidth")

  kern = match.arg(kern)
  if(kern == "triangle"){

    res <- pmax(1 - abs(x - t) / h, 0)
  }

  return(res)
}
