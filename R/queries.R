#' @title Whether a 'ratioOfQsprays' is polynomial
#' @description Checks whether a \code{ratioOfQsprays} actually is polynomial,
#'   that is, whether its denominator is a constant \code{qspray} polynomial
#'   (and then it should be equal to one).
#'
#' @param roq a \code{ratioOfQsprays} object
#'
#' @return A Boolean value.
#' @export
#' @importFrom qspray isConstantQspray
#'
#' @examples
#' x <- qlone(1)
#' y <- qlone(2)
#' roq <- (x^2 - y^2) / (x - y)
#' isPolynomial(roq)
#' roq == x + y
isPolynomial <- function(roq) {
  isConstantQspray(roq@denominator)
}
