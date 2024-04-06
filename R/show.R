#' @title Print a 'ratioOfQsprays' object
#' @description Print a \code{ratioOfQsprays} object given a function to print
#'   a \code{qspray} object
#'
#' @param showQspray a function which prints a \code{qspray} object, which will
#'   be applied to the numerator and the denominator
#' @param quotientBar a string representing the quotient bar between the
#'   numerator and the denominator, including surrounding spaces,
#'   e.g \code{" / "}
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#'
#' @seealso \code{\link{showRatioOfQspraysCanonical}}
showRatioOfQsprays <- function(showQspray, quotientBar = "  %//%  ") {
  function(roq) {
    if(isQone(roq@denominator)) {
      sprintf(
        "[%s]", showQspray(roq@numerator)
      )
    } else {
      sprintf(
        "[ %s ]%s[ %s ]",
        showQspray(roq@numerator),
        quotientBar,
        showQspray(roq@denominator)
      )
    }
  }
}

#' @title Print a 'ratioOfQsprays'
#' @description Print a \code{ratioOfQsprays} object given a string to denote
#'   the unindexed variables.
#'
#' @param var a string, usually a letter, to denote the unindexed variables
#' @param quotientBar a string representing the quotient bar between the
#'   numerator and the denominator, including surrounding spaces,
#'   e.g \code{" / "}
#' @param ... arguments other than \code{quotientBar} passed to
#'   \code{\link{showRatioOfQsprays}} (currently there's no such argument)
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#' @importFrom qspray showQsprayCanonical
#'
#' @note The \code{show} method for \code{ratioOfQsprays} objects uses
#'   \code{showRatioOfQspraysCanonical("x", quotientBar = "  \%//\%  ")}
#'   by default. But this can be controlled as follows. If a
#'   \code{ratioOfQsprays} object has an attribute \code{"x"}, then the value
#'   of this attribute will replace \code{"x"} in the \code{show} output.
#'   It is also possible to control the \code{quotientBar} argument by
#'   assigning a \code{"quotientBar"} attribute to the \code{ratioOfQsprays}
#'   object to be printed.
#'
#' @examples
#' roq <- rRatioOfQsprays
#' showRatioOfQspraysCanonical("X", " / ")(roq)
showRatioOfQspraysCanonical <- function(var, quotientBar = "  %//%  ", ...) {
  showRatioOfQsprays(showQsprayCanonical(var), quotientBar = quotientBar, ...)
}
