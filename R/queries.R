#' @include ratioOfQsprays.R
NULL

setGeneric("numberOfVariables")
setGeneric("isConstant")
setGeneric("isUnivariate")

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

#' @name isUnivariate
#' @aliases isUnivariate,ratioOfQsprays-method
#' @docType methods
#' @importFrom qspray isUnivariate
#' @title Whether a 'ratioOfQsprays' is univariate
#' @description Checks whether a \code{ratioOfQsprays} object defines a
#'   univariate fraction of polynomials.
#'
#' @param x a \code{ratioOfQsprays} object
#'
#' @return A Boolean value.
#' @export
setMethod(
  "isUnivariate", "ratioOfQsprays",
  function(x) {
    numberOfVariables(x) %in% c(0L, 1L)
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

#' @title Get the numerator of a 'ratioOfQsprays'
#' @description Get the numerator of a \code{ratioOfQsprays} object,
#'   preserving the show options.
#'
#' @param roq a \code{ratioOfQsprays} object
#'
#' @return A \code{qspray} object.
#' @export
getNumerator <- function(roq) {
  qspray <- roq@numerator
  passShowAttributes(roq, qspray)
}

#' @title Get the denominator of a 'ratioOfQsprays'
#' @description Get the denominator of a \code{ratioOfQsprays} object,
#'   preserving the show options.
#'
#' @param roq a \code{ratioOfQsprays} object
#'
#' @return A \code{qspray} object.
#' @export
getDenominator <- function(roq) {
  qspray <- roq@denominator
  passShowAttributes(roq, qspray)
}
