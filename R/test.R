#' Add two numbers
#'
#' Calculate the sum of two numbers
#'
#' @param a first number
#' @param b second number
#'
#' @return the sum of two numbers
#' @export
#'
#' @examples
#' test_fun(1, 3)
test_fun <- function(a, b){

  a + b
}


#' RDHonest result for Lee dataset
#'
#' Return point estimate for RDHonest procedure on Lee dataset.
#'
#' This function was made for a test purpose for package import.
#'
#' @return point estimate for RDHonest procedure on Lee dataset
#' @export
#'
#' @examples
#' RDH_lee()
RDH_lee <- function(){

  res <- RDHonest::RDHonest(voteshare ~ margin, data = RDHonest::lee08, kern = "uniform", M = 0.1, h = 10, sclass = "T")
  return(res$estimate)
}
