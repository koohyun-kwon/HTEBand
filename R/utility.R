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
