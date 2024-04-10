#' @title Print a 'ratioOfQsprays' object
#' @description Prints a \code{ratioOfQsprays} object given a function to print
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
#' @seealso \code{\link{showRatioOfQspraysX1X2X3}},
#'   \code{\link{showRatioOfQspraysXYZ}}.
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
#'   the non-indexed variables.
#'
#' @param var a string, usually a letter, to denote the non-indexed variables
#' @param quotientBar a string representing the quotient bar between the
#'   numerator and the denominator, including surrounding spaces,
#'   e.g \code{" / "}
#' @param ... arguments other than \code{quotientBar} passed to
#'   \code{\link{showRatioOfQsprays}} (currently there's no such argument)
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#' @importFrom qspray showQsprayX1X2X3
#'
#' @examples
#' ( roq <- rRatioOfQsprays() )
#' showRatioOfQspraysX1X2X3("X", " / ")(roq)
showRatioOfQspraysX1X2X3 <- function(var, quotientBar = "  %//%  ", ...) {
  showRatioOfQsprays(showQsprayX1X2X3(var), quotientBar = quotientBar, ...)
}

#' @title Print a 'ratioOfQsprays'
#' @description Print a \code{ratioOfQsprays} object given some letters to
#'   denote the variables, by printing monomials like \code{"x^2yz"}.
#'
#' @param letters a vector of strings, usually some letters such as \code{"x"}
#'   and \code{"y"}, to denote the variables
#' @param quotientBar a string representing the quotient bar between the
#'   numerator and the denominator, including surrounding spaces,
#'   e.g \code{" / "}
#' @param ... arguments other than \code{quotientBar} passed to
#'   \code{\link{showRatioOfQsprays}} (currently there's no such argument)
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#' @importFrom qspray showQsprayXYZ
#'
#' @examples
#' ( roq <- rRatioOfQsprays() )
#' showRatioOfQspraysXYZ(c("X", "Y", "Z"), " / ")(roq)
showRatioOfQspraysXYZ <- function(letters, quotientBar = "  %//%  ", ...) {
  showRatioOfQsprays(showQsprayXYZ(letters), quotientBar = quotientBar, ...)
}

#' @title Set show option to a 'qspray' object
#' @description Set show option to a \code{qspray} object
#'
#' @param x a \code{qspray} object
#' @param which which option to set; this can be \code{"x"},
#'   \code{"quotientBar"}, \code{"showQspray"}, or \code{"showRatioOfQsprays"}
#' @param value the value of the option to be set
#'
#' @return This returns the updated \code{qspray}.
#' @export
#'
#' @examples
#' ( roq <- rRatioOfQsprays() )
#' showRatioOfQspraysOption(roq, "quotientBar") <- " / "
#' roq
#' showRatioOfQspraysOption(roq, "x") <- "a"
#' roq
`showRatioOfQspraysOption<-` <- function(x, which, value) {
  which <-
    match.arg(which, c("x", "quotientBar", "showQspray", "showRatioOfQsprays"))
  showOpts <- attr(x, "showOpts") %||% TRUE
  attr(showOpts, which) <- value
  if(which != "showRatioOfQsprays") {
    if(which == "x") {
      sQ <- showQsprayX1X2X3(value)
    } else if(which == "quotientBar") {
      trivariate <- numberOfVariables(x) <= 3L
      if(trivariate) {
        sQ <- showQsprayXYZ()
      } else {
        sQ <- showQsprayX1X2X3(attr(showOpts, "x") %||% "x")
      }
    } else if(which == "showQspray") {
      sQ <- value
    }
    sROQ <- showRatioOfQsprays(
      showQspray = sQ,
      quotientBar = attr(showOpts, "quotientBar") %||% "  %//%  "
    )
  } else {
    sQ   <- NULL
    sROQ <- value
  }
  attr(showOpts, "showQspray") <- sQ
  attr(showOpts, "showRatioOfQsprays") <- sROQ
  attr(x, "showOpts") <- showOpts
  x
}

getShowRatioOfQsprays <- function(roq) {
  showOpts <- attr(roq, "showOpts")
  attr(showOpts, "showRatioOfQsprays") %||%
    attr(attr(showOpts, "showSymbolicQspray"), "showRatioOfQsprays") %||%
    showRatioOfQsprays(
      showQspray = attr(showOpts, "showQspray") %||%
        showQsprayX1X2X3(attr(showOpts, "x") %||% "x"),
      quotientBar = attr(showOpts, "quotientBar") %||% "  %//%  "
    )
}
