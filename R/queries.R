setGeneric("numberOfVariables")
setGeneric("isConstant")

#' @name numberOfVariables
#' @aliases numberOfVariables,ratioOfQsprays-method
#' @docType methods
#' @importFrom qspray numberOfVariables
#' @title Number of variables in a 'ratioOfQsprays'
#' @description Number of variables involved in a \code{ratioOfQsprays} object.
#'
#' @param x a \code{ratioOfQsprays} object
#'
#' @return An integer.
#' @export
#' @note The number of variables in the \code{qspray} object \code{qlone(d)}
#'   is \code{d}, not \code{1}.
setMethod(
  "numberOfVariables", "ratioOfQsprays",
  function(x) {
    max(numberOfVariables(x@numerator), numberOfVariables(x@denominator))
  }
)

#' @name isConstant
#' @aliases isConstant,ratioOfQsprays-method
#' @docType methods
#' @importFrom qspray isConstant
#' @title Whether a 'ratioOfQsprays' is constant
#' @description Checks whether a \code{ratioOfQsprays} object defines a constant
#'   fraction of polynomials.
#'
#' @param x a \code{ratioOfQsprays} object
#'
#' @return A Boolean value.
#' @export
setMethod(
  "isConstant", "ratioOfQsprays",
  function(x) {
    numberOfVariables(x) == 0L
  }
)

#' @title Whether a 'ratioOfQsprays' is polynomial
#' @description Checks whether a \code{ratioOfQsprays} actually is polynomial,
#'   that is, whether its denominator is a constant \code{qspray} polynomial
#'   (and then it should be equal to one).
#'
#' @param roq a \code{ratioOfQsprays} object
#'
#' @return A Boolean value.
#' @export
#' @importFrom qspray isConstant
#'
#' @examples
#' x <- qlone(1)
#' y <- qlone(2)
#' roq <- (x^2 - y^2) / (x - y)
#' isPolynomial(roq)
#' roq == x + y
isPolynomial <- function(roq) {
  isConstant(roq@denominator)
}
