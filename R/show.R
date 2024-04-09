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
#' @importFrom qspray showQsprayCanonical showQsprayUnivariate
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
#' roq <- rRatioOfQsprays()
#' showRatioOfQspraysCanonical("X", " / ")(roq)
showRatioOfQspraysCanonical <- function(var, quotientBar = "  %//%  ", ...) {
  showRatioOfQsprays(showQsprayCanonical(var), quotientBar = quotientBar, ...)
}

#' Title
#'
#' @param var xx
#' @param quotientBar xx
#' @param ... xx
#'
#' @return xx
#' @export
showRatioOfQspraysUnivariate <- function(var, quotientBar = "  %//%  ", ...) {
  showRatioOfQsprays(showQsprayUnivariate(var), quotientBar = quotientBar, ...)
}

#' @title Set show options to a 'ratioOfQsprays' object
#' @description Set some attributes to a \code{ratioOfQsprays} object
#'   to control the way it is displayed. See the note in the
#'   documentation of \code{\link{showRatioOfQspraysCanonical}} for details.
#'
#' @param roq a \code{ratioOfQsprays} object
#' @param x value for the \code{"x"} attribute
#' @param quotientBar value for the \code{"quotientBar"} attribute
#'
#' @return The input \code{ratioOfQsprays} object with new attributes.
#' @export
#'
#' @examples
#' ( roq <- rRatioOfQsprays() )
#' withAttributes(roq, x = "a", quotientBar = " / ")
withAttributes <- function(
    roq, x = "x", quotientBar = "  %//%  "
) {
  attr(roq, "x") <- x
  attr(roq, "quotientBar") <- quotientBar
  roq
}

#' @title Set show option to a 'qspray' object
#' @description Set show option to a \code{qspray} object
#'
#' @param x a \code{qspray} object
#' @param which which option to set; this can be \code{"x"},
#'   \code{"showMonomial"}, or \code{"showQspray"}
#' @param value the value of the option
#'
#' @return This returns the updated \code{qspray}.
#' @export
#'
#' @examples
#' roq <- rRatioOfQsprays()
#' showRatioOfQspraysOption(roq, "x") <- "a"
#' showRatioOfQspraysOption(roq, "quotientBar") <- " / "
#' roq
`showRatioOfQspraysOption<-` <- function(x, which, value) {
  which <-
    match.arg(which, c("x", "quotientBar", "showQspray", "showRatioOfQsprays"))
  showOpts <- attr(x, "showOpts") %||% TRUE
  attr(showOpts, which) <- value
  univariate <- numberOfVariables2(x) == 1L
  sROQ <- if(univariate) {
    showRatioOfQspraysUnivariate
  } else {
    showRatioOfQspraysCanonical
  }
  if(which == "x") {
    attr(showOpts, "showRatioOfQsprays") <-
      sROQ(
        var = value,
        quotientBar = attr(showOpts, "quotientBar") %||% "  %//%  "
      )
  } else if(which == "quotientBar") {
    attr(showOpts, "showRatioOfQsprays") <-
      sROQ(
        var = attr(showOpts, "x") %||% "x",
        quotientBar = value
      )
  } else if(which == "showQspray") {
    attr(showOpts, "showRatioOfQsprays") <- showRatioOfQsprays(
      showQspray = value,
      quotientBar = attr(showOpts, "quotientBar") %||% "  %//%  "
    )
  }
  attr(x, "showOpts") <- showOpts
  x
}

getShowRatioOfQsprays <- function(roq) {
  showOpts <- attr(roq, "showOpts")
  attr(showOpts, "showRatioOfQsprays") %||%
    attr(attr(showOpts, "showSymbolicQspray"), "showRatioOfQsprays") %||%
    showRatioOfQsprays(
      showQspray = attr(showOpts, "showQspray") %||%
        showQsprayCanonical(attr(showOpts, "a") %||% "a"),
      quotientBar = attr(showOpts, "quotientBar") %||% " %//% "
    )
}
