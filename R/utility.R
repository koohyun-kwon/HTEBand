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
