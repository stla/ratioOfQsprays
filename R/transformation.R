#' @title Partial derivative
#' @description Partial derivative of a \code{ratioOfQsprays}.
#'
#' @param roq object of class \code{ratioOfQsprays}
#' @param i integer, the dimension to differentiate with respect to
#' @param derivative integer, how many times to differentiate
#'
#' @return A \code{ratioOfQsprays} object.
#' @importFrom qspray derivQspray
#' @export
#'
#' @examples
#' library(ratioOfQsprays)
#' x <- qlone(1)
#' y <- qlone(2)
#' roq <- (2*x  + 3*x*y) / (x^2 + y^2)
#' derivRatioOfQsprays(roq, 1)
derivRatioOfQsprays <- function(roq, i, derivative = 1) {
  stopifnot(inherits(roq, "ratioOfQsprays"))
  stopifnot(isNonnegativeInteger(i))
  stopifnot(isPositiveInteger(derivative))
  f  <- roq@numerator
  g  <- roq@denominator
  fp <- derivQspray(f, i)
  gp <- derivQspray(g, i)
  (fp*g - f*gp) / g^2
}

#' @title Partial differentiation
#' @description Partial differentiation of a ratioOfQsprays polynomial.
#'
#' @param roq object of class \code{ratioOfQsprays}
#' @param orders integer vector, the orders of the differentiation
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#'
#' @examples
#' library(ratioOfQsprays)
#' x <- qlone(1)
#' y <- qlone(2)
#' roq <- (x + 2*y  + 3*x*y) / (x + 1)
#' dRatioOfQsprays(roq, c(1, 1))
#' derivRatioOfQsprays(derivRatioOfQsprays(roq, 1), 2)
dRatioOfQsprays <- function(roq, orders) {
  stopifnot(inherits(roq, "ratioOfQsprays"))
  for(i in seq_along(orders)) {
    stopifnot(isPositiveInteger(orders[i]))
  }
  orders <- removeTrailingZeros(orders)
  if(length(orders) > numberOfVariables2(roq)) {
    return(as.ratioOfQsprays(0L))
  }
  ns <- do.call(c, lapply(seq_along(as.integer(orders)), function(i) {
    rep(i, orders[i])
  }))
  f <- function(r, i) {
    if(i != 0L) {
      derivRatioOfQsprays(r, i)
    } else {
      r
    }
  }
  Reduce(f, ns, init = roq)
}

#' @title Permute variables
#' @description Permute the variables of a \code{ratioOfQsprays} polynomial.
#'
#' @param roq a \code{ratioOfQsprays} object
#' @param permutation a permutation
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#' @importFrom qspray permuteVariables
#'
#' @examples
#' library(ratioOfQsprays)
#' f <- function(x, y, z) {
#'   (x^2 + 5*y + z - 1) / (x + 1)
#' }
#' x <- qlone(1)
#' y <- qlone(2)
#' z <- qlone(3)
#' R <- f(x, y, z)
#' permutation <- c(3, 1, 2)
#' S <- permuteVariables2(R, permutation)
#' S == f(z, x, y) # should be TRUE
permuteVariables2 <- function(roq, permutation) {
  permuteVariables(roq@numerator, permutation) /
    permuteVariables(roq@denominator, permutation)
}

#' @title Swap variables
#' @description Swap two variables of a \code{ratioOfQsprays}.
#'
#' @param roq a \code{ratioOfQsprays} object
#' @param i,j indices of the variables to be swapped
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#' @importFrom qspray swapVariables
#'
#' @examples
#' library(ratioOfQsprays)
#' f <- function(x, y, z) {
#'   (x^2 + 5*y + z - 1) / (x + 1)
#' }
#' x <- qlone(1)
#' y <- qlone(2)
#' z <- qlone(3)
#' R <- f(x, y, z)
#' S <- swapVariables2(R, 2, 3)
#' S == f(x, z, y) # should be TRUE
swapVariables2 <- function(roq, i, j) {
  swapVariables(roq@numerator, i, j) /
    swapVariables(roq@denominator, i, j)
  # new(
  #   "ratioOfQsprays",
  #   numerator   = swapVariables(roq@numerator, i, j),
  #   denominator = swapVariables(roq@denominator, i, j)
  # )
}
