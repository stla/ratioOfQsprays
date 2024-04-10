#' @title Print a 'ratioOfQsprays' object
#' @description Prints a \code{ratioOfQsprays} object given a function to print
#'   a \code{qspray} object
#'
#' @param showQspray a function which prints a \code{qspray} object, which will
#'   be applied to the numerator and the denominator
#' @param quotientBar a string representing the quotient bar between the
#'   numerator and the denominator, including surrounding spaces,
#'   e.g \code{" / "}
#' @param lbracket,rbracket used to enclose the numerator and the denominator
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#'
#' @seealso \code{\link{showRatioOfQspraysX1X2X3}},
#'   \code{\link{showRatioOfQspraysXYZ}},
#'   \code{\link{showRatioOfQspraysOption<-}}.
#'
#' @examples
#' set.seed(666)
#' ( roq <- rRatioOfQsprays() )
#' f <- showRatioOfQsprays(showQsprayX1X2X3("a"), " / ", "[[[ ", " ]]]")
#' f(roq)
showRatioOfQsprays <- function(
  showQspray, quotientBar = "  %//%  ", lbracket = "[ ", rbracket = " ]"
) {
  function(roq) {
    enclose <- function(qspray) {
      sprintf(
        "%s%s%s", lbracket, showQspray(qspray), rbracket
      )
    }
    if(isQone(roq@denominator)) {
      enclose(roq@numerator)
    } else {
      sprintf(
        "%s%s%s",
        enclose(roq@numerator),
        quotientBar,
        enclose(roq@denominator)
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
#'   \code{\link{showRatioOfQsprays}}
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#' @importFrom qspray showQsprayX1X2X3
#'
#' @seealso \code{\link{showRatioOfQspraysXYZ}},
#'   \code{\link{showRatioOfQspraysOption<-}}.
#'
#' @examples
#' set.seed(666)
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
#'   \code{\link{showRatioOfQsprays}}
#'
#' @return A function which takes as argument a \code{ratioOfQsprays} object
#'   and which prints it.
#' @export
#' @importFrom qspray showQsprayXYZ
#'
#' @seealso \code{\link{showRatioOfQspraysX1X2X3}},
#'   \code{\link{showRatioOfQspraysOption<-}}.
#'
#' @examples
#' set.seed(666)
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
#' @importFrom qspray showQspray showMonomialXYZ showMonomialX1X2X3
#'
#' @examples
#' set.seed(666)
#' ( roq <- rRatioOfQsprays() )
#' showRatioOfQspraysOption(roq, "quotientBar") <- " / "
#' roq
#' showRatioOfQspraysOption(roq, "x") <- "a"
#' roq
#' showRatioOfQspraysOption(roq, "showQspray") <- showQsprayXYZ()
#' roq
`showRatioOfQspraysOption<-` <- function(x, which, value) {
  which <-
    match.arg(which, c("x", "quotientBar", "showQspray", "showRatioOfQsprays"))
  showOpts <- attr(x, "showOpts") %||% TRUE
  attr(showOpts, which) <- value
  if(which == "inheritable") {
    attr(x, "showOpts") <- showOpts
    return(x)
  }
  if(which != "showRatioOfQsprays") {
    if(which == "x") {
      univariate <- isUnivariate(x)
      if(univariate) {
        sM <- showMonomialXYZ(letters = value)
      } else {
        sM <- showMonomialX1X2X3(x = value)
        attr(showOpts, "inheritable") <- TRUE
      }
      attr(showOpts, "showMonomial") <- sM
      sQ <- showQspray(sM)
    } else if(which == "quotientBar") {
      trivariate <- numberOfVariables(x) <= 3L
      if(trivariate) {
        sM <- showMonomialXYZ()
      } else {
        sM <- showMonomialX1X2X3(attr(showOpts, "x") %||% "x")
        attr(showOpts, "inheritable") <- TRUE
      }
      sQ <- showQspray(sM)
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

setDefaultShowRatioOfQspraysOption <- function(roq) {
  showRatioOfQspraysOption(roq, "quotientBar") <- "  %//%  "
  invisible(roq)
}

getShowRatioOfQsprays <- function(roq) {
  showOpts <- attr(roq, "showOpts")
  sROQ <- attr(showOpts, "showRatioOfQsprays")
  if(is.null(sROQ)) {
    # it's possible that showOpts has a "showQspray" attribute, from a call
    # to as.RatioOfQsprays
    sQ <- attr(showOpts, "showQspray")
    if(is.null(sQ)) {
      roq <- setDefaultShowRatioOfQspraysOption(roq)
      sROQ <- attr(attr(roq, "showOpts"), "showRatioOfQsprays")
    } else {
      showRatioOfQspraysOption(roq, "showQspray") <- sQ
    }
  }
  attr(attr(roq, "showOpts"), "showRatioOfQsprays")
}
