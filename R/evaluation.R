#' @title Evaluate a 'ratioOfQsprays' object
#' @description Evaluation of the fraction of multivariate polynomials
#'  represented by a \code{ratioOfQsprays} object.
#'
#' @param roq a \code{ratioOfQsprays} object
#' @param values_re vector of the real parts of the values; each element of
#'   \code{as.character(values_re)} must be a quoted integer or a quoted fraction
#' @param values_im vector of the imaginary parts of the values; each element of
#'   \code{as.character(values_im)} must be a quoted integer or a quoted fraction
#'
#' @return A \code{bigq} number if \code{values_im=NULL}, a pair of \code{bigq}
#'   numbers otherwise: the real part and the imaginary part of the result.
#' @export
#' @importFrom qspray evalQspray
#'
#' @examples
#' x <- qlone(1); y <- qlone(2)
#' roq <- 2*x / (x^2 + 3*y^2)
#' evalRatioOfQsprays(roq, c("2", "5/2", "99999")) # "99999" will be ignored
evalRatioOfQsprays <- function(roq, values_re, values_im = NULL) {
  num <- evalQspray(roq@numerator, values_re, values_im)
  den <- evalQspray(roq@denominator, values_re, values_im)
  if(is.null(values_im)) {
    num / den
  } else {
    a <- num[1L]
    b <- num[2L]
    c <- den[1L]
    d <- den[2L]
    den  <- c*c + d*d
    numRe <- a*c + b*d
    numIm <- b*c - a*d
    c(numRe/den, numIm/den)
  }
}

#' @title Substitute the variables of a 'ratioOfQsprays' with 
#'   some 'ratioOfQsprays' 
#' @description Substitute the variables of a \code{ratioOfQsprays} fraction 
#'   of polynomials with some \code{ratioOfQsprays} objects, or some objects 
#'   coercible to some \code{ratioOfQsprays}.
#'
#' @param roq a \code{ratioOfQsprays} object
#' @param listOfRoqs the \code{ratioOfQsprays} to be substituted to the 
#'   variables of \code{roq}; this must be a list of some objects coercible 
#'   to some \code{ratioOfQsprays} objects (by an application of the function
#'   \code{\link{as.ratioOfQsprays}}), and the length of this list must be 
#'   equal to the number of variables of \code{roq}
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#' @examples
#' x <- qlone(1)
#' y <- qlone(2)
#' z <- qlone(3)
#' p <- x + y
#' q <- x - y
#' rOQ <- p / q
#' rOQ1 <- x / y
#' rOQ2 <- z / x
#' substituteSomeRatioOfQsprays(rOQ, list(rOQ1, rOQ2))
substituteSomeRatioOfQsprays <- function(roq, listOfRoqs) {
  stopifnot(is.list(listOfRoqs))
  n <- numberOfVariables(roq)
  if(n != length(listOfRoqs)) {
    stop(
      "The length of the `listOfRoqs` list is not correct."
    ) 
  }
  listOfRoqs <- lapply(listOfRoqs, as.ratioOfQsprays)
  qspray1 <- roq@numerator
  powers1 <- qspray1@powers
  coeffs1 <- as.bigq(qspray1@coeffs)
  qspray2 <- roq@denominator
  powers2 <- qspray2@powers
  coeffs2 <- as.bigq(qspray2@coeffs)
  zeroRatioOfQsprays <- as.ratioOfQsprays(0L)
  num <- zeroRatioOfQsprays
  for(i in seq_along(powers1)) {
    exponents <- powers1[[i]]
    term <- 1L
    for(j in seq_along(exponents)) {
      term <- term * listOfRoqs[[j]]^exponents[j]
    }
    num <- num + coeffs1[i] * term
  }
  den <- zeroRatioOfQsprays
  for(i in seq_along(powers2)) {
    exponents <- powers2[[i]]
    term <- 1L
    for(j in seq_along(exponents)) {
      term <- term * listOfRoqs[[j]]^exponents[j]
    }
    den <- den + coeffs2[i] * term
  }
  num / den
}


#' @title Partial evaluation of a 'ratioOfQsprays' fraction of polynomials
#' @description Substitute some values to a subset of the variables of a
#'   \code{ratioOfQsprays} fraction of polynomials.
#'
#' @param roq a \code{ratioOfQsprays} object
#' @param values the values to be substituted; this must be a vector whose
#'   length equals the number of variables of \code{roq}, and whose each
#'   entry is either \code{NA} for non-substitution or a "scalar" \code{x}
#'   such that \code{as.character(x)} is a quoted integer or a quoted fraction,
#'   e.g. a \code{bigq} number
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#' @importFrom qspray substituteQspray numberOfVariables
#' @importFrom utils head
#'
#' @examples
#' library(ratioOfQsprays)
#' x <- qlone(1)
#' y <- qlone(2)
#' z <- qlone(3)
#' roq <- (x^2 + y^2 + x*y*z - 1) / (x + 1)
#' substituteRatioOfQsprays(roq, c("2", NA, "3/2"))
substituteRatioOfQsprays <- function(roq, values) {
  n1 <- numberOfVariables(roq@numerator)
  n2 <- numberOfVariables(roq@denominator)
  substituteQspray(roq@numerator, head(values, n1)) /
    substituteQspray(roq@denominator, head(values, n2))
}

#' @title Ratio of multivariate polynomials as function
#' @description Coerces a \code{ratioOfQsprays} fraction of polynomials 
#'   to a function.
#'
#' @param x object of class \code{ratioOfQsprays}
#' @param N Boolean, whether the function must numerically approximate
#'   the result
#' @param ... ignored
#'
#' @return A function having the same variables as \code{x}. If
#'   \code{N=FALSE}, this function returns a string. If \code{N=TRUE}, it 
#'   returns a number if the result does not contain any variable, otherwise 
#'   it returns a R expression.
#' @export
#' @importFrom Ryacas yac_str as_r
#' @importFrom methods formalArgs
#'
#' @examples
#' library(ratioOfQsprays)
#' x <- qlone(1); y <- qlone(2)
#' roq <- (x^2/2 + y^2 + x*y - 1) / (x + 1)
#' f <- as.function(roq)
#' g <- as.function(roq, N = TRUE)
#' f(2, "3/7")
#' g(2, "3/7")
#' f("x", "y")
#' g("x", "y")
#' # the evaluation is performed by (R)yacas and complex numbers are
#' # allowed; the imaginary unit is denoted by \code{I}:
#' f("2 + 2*I", "Sqrt(2)")
#' g("2 + 2*I", "Sqrt(2)")
as.function.ratioOfQsprays <- function(x, N = FALSE, ...) {
  fnum <- as.function(x@numerator, N = FALSE)
  fden <- as.function(x@denominator, N = FALSE)
  formalsNum <- formals(fnum)
  formalsDen <- formals(fden)
  if(length(formalsNum) > length(formalsDen)) {
    formals(fden) <- formalsNum
  } else {
    formals(fnum) <- formalsDen
  }
  vars <- formalArgs(fnum)
  if(N) {
    f <- function() {
      do.call(function(...) {
        as_r(yac_str(
          sprintf(
            "N((%s)/(%s))",
            as.character(fnum(...)), as.character(fden(...))
          )
        ))
      }, lapply(vars, function(xi) {
        eval(parse(text = xi))
      }))
    }
  } else {
    f <- function() {
      do.call(function(...) {
        yac_str(
          sprintf("(%s)/(%s)", as.character(fnum(...)), as.character(fden(...)))
        )
      }, lapply(vars, function(xi) {
        eval(parse(text = xi))
      }))
    }
  }
  formals(f) <- formals(fnum)
  f
}
