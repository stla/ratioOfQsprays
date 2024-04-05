#' @useDynLib ratioOfQsprays, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom qspray qone as.qspray qsprayDivision isQone isQzero isConstantQspray getConstantTerm
#' @importFrom resultant gcd
#' @importFrom methods setMethod setClass new show
#' @importFrom gmp as.bigq
#' @include ratioOfQsprays.R
NULL

setClass(
  "ratioOfQsprays",
  slots = c(numerator = "qspray", denominator = "qspray")
)

#' @importFrom utils capture.output
#' @noRd
showRatioOfQsprays <- function(roq) {
  # if(roq@numerator == qzero()) {
  #   return("0")
  # }
  if(isQone(roq@denominator)) {
    sprintf(
      "[%s]",
      trimws(capture.output(show(roq@numerator)),   which = "right")
    )
  } else {
    sprintf(
      "[%s] / [%s]",
      trimws(capture.output(show(roq@numerator)),   which = "right"),
      trimws(capture.output(show(roq@denominator)), which = "right")
    )
  }
}

setMethod(
  "show", "ratioOfQsprays",
  function(object) {
    cat(showRatioOfQsprays(object), "\n")
  }
)

identifyQspray <- function(qspray) {
  new("ratioOfQsprays", numerator = qspray, denominator = qone())
}

setGeneric(
  "as.ratioOfQsprays", function(x) {
    NULL
  }
)

#' @name as.ratioOfQsprays
#' @aliases as.ratioOfQsprays,character-method as.ratioOfQsprays,ratioOfQsprays-method as.ratioOfQsprays,qspray-method as.ratioOfQsprays,numeric-method as.ratioOfQsprays,bigz-method as.ratioOfQsprays,bigq-method
#' @exportMethod as.ratioOfQsprays
#' @docType methods
#' @title Coercion to a 'ratioOfQsprays' object
#'
#' @param x a \code{ratioOfQsprays} object, a \code{qspray} object, or an
#'   object yielding a quoted integer or a quoted fraction after an application
#'   of \code{as.character}
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#'
#' @examples
#' library(qspray)
#' as.ratioOfQsprays(2)
#' as.ratioOfQsprays("1/3")
#' as.ratioOfQsprays(5*qlone(1) + qlone(2)^2)
setMethod(
  "as.ratioOfQsprays", "character",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "ratioOfQsprays",
  function(x) {
    x
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "qspray",
  function(x) {
    identifyQspray(x)
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "numeric",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "bigz",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "bigq",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @name ratioOfQsprays-unary
#' @title Unary operators for ratioOfQsprays objects
#' @description Unary operators for ratioOfQsprays objects.
#' @aliases +,ratioOfQsprays,missing-method -,ratioOfQsprays,missing-method
#' @param e1 object of class \code{ratioOfQsprays}
#' @param e2 nothing
#' @return A \code{ratioOfQsprays} object.
setMethod(
  "+",
  signature(e1 = "ratioOfQsprays", e2 = "missing"),
  function(e1, e2) e1
)
#' @rdname ratioOfQsprays-unary
setMethod(
  "-",
  signature(e1 = "ratioOfQsprays", e2 = "missing"),
  function(e1, e2) {
    new(
      "ratioOfQsprays",
      powers = e1@powers, coeffs = as.character(-as.bigq(e1@coeffs))
    )
  }
)

simplifyRatioOfQsprays <- function(roq) {
  num <- roq@numerator
  den <- roq@denominator
  g <- gcd(num, den)
  num <- qsprayDivision(num, g)[["Q"]]
  den <- qsprayDivision(den, g)[["Q"]]
  if(isConstantQspray(den)) {
    k <- 1L / getConstantTerm(den)
    num <- k * num
    den <- k * den
  }
  new(
    "ratioOfQsprays",
    numerator   = num,
    denominator = den
  )
}

ratioOfQsprays_arith_ratioOfQsprays <- function(e1, e2) {
  num1 <- e1@numerator
  den1 <- e1@denominator
  num2 <- e2@numerator
  den2 <- e2@denominator
  switch(
    .Generic,
    "+" = {
      x <- ROQaddition(
        list("powers" = num1@powers, "coeffs" = num1@coeffs),
        list("powers" = den1@powers, "coeffs" = den1@coeffs),
        list("powers" = num2@powers, "coeffs" = num2@coeffs),
        list("powers" = den2@powers, "coeffs" = den2@coeffs)
      )
      simplifyRatioOfQsprays(
        new(
          "ratioOfQsprays",
          numerator   = qspray_from_list(x[["numerator"]]),
          denominator = qspray_from_list(x[["denominator"]])
        )
      )
    },
    "-" = {
      x <- ROQsubtraction(
        list("powers" = num1@powers, "coeffs" = num1@coeffs),
        list("powers" = den1@powers, "coeffs" = den1@coeffs),
        list("powers" = num2@powers, "coeffs" = num2@coeffs),
        list("powers" = den2@powers, "coeffs" = den2@coeffs)
      )
      simplifyRatioOfQsprays(
        new(
          "ratioOfQsprays",
          numerator   = qspray_from_list(x[["numerator"]]),
          denominator = qspray_from_list(x[["denominator"]])
        )
      )
    },
    "*" = {
      x <- ROQmultiplication(
        list("powers" = num1@powers, "coeffs" = num1@coeffs),
        list("powers" = den1@powers, "coeffs" = den1@coeffs),
        list("powers" = num2@powers, "coeffs" = num2@coeffs),
        list("powers" = den2@powers, "coeffs" = den2@coeffs)
      )
      simplifyRatioOfQsprays(
        new(
          "ratioOfQsprays",
          numerator   = qspray_from_list(x[["numerator"]]),
          denominator = qspray_from_list(x[["denominator"]])
        )
      )
    },
    "/" = {
      x <- ROQdivision(
        list("powers" = num1@powers, "coeffs" = num1@coeffs),
        list("powers" = den1@powers, "coeffs" = den1@coeffs),
        list("powers" = num2@powers, "coeffs" = num2@coeffs),
        list("powers" = den2@powers, "coeffs" = den2@coeffs)
      )
      simplifyRatioOfQsprays(
        new(
          "ratioOfQsprays",
          numerator   = qspray_from_list(x[["numerator"]]),
          denominator = qspray_from_list(x[["denominator"]])
        )
      )
    },
    stop(gettextf(
      "Binary operator %s not defined for ratioOfQsprays objects.",
      dQuote(.Generic)
    ))
  )
}
ratioOfQsprays_arith_qspray <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = e1 * as.ratioOfQsprays(e2),
    "/" = e1 / as.ratioOfQsprays(e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
qspray_arith_ratioOfQsprays <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.ratioOfQsprays(e1) + e2,
    "-" = as.ratioOfQsprays(e1) - e2,
    "*" = as.ratioOfQsprays(e1) * e2,
    "/" = as.ratioOfQsprays(e1) / e2,
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
ratioOfQsprays_arith_character <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = new(
      "ratioOfQsprays",
      numerator   = e1@numerator * e2,
      denominator = e1@denominator
    ),
    "/" = new(
      "ratioOfQsprays",
      numerator   = e1@numerator / e2,
      denominator = e1@denominator
    ),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
ratioOfQspraysPower <- function(ratioOfQsprays, n) {
  stopifnot(isInteger(n))
  numerator   <- ratioOfQsprays@numerator
  denominator <- ratioOfQsprays@denominator
  roqAsList <- ROQpower(
    list("powers" = numerator@powers,   "coeffs" = numerator@coeffs),
    list("powers" = denominator@powers, "coeffs" = denominator@coeffs),
    n
  )
  simplifyRatioOfQsprays(
    new(
      "ratioOfQsprays",
      numerator   = qspray_from_list(roqAsList[["numerator"]]),
      denominator = qspray_from_list(roqAsList[["denominator"]])
    )
  )
}
ratioOfQsprays_arith_gmp <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = new(
      "ratioOfQsprays",
      numerator = e1@numerator * e2,
      denominator = e1@denominator
    ),
    "/" = new(
      "ratioOfQsprays",
      numerator   = e1@numerator / e2,
      denominator = e1@denominator
    ),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
ratioOfQsprays_arith_numeric <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = new(
      "ratioOfQsprays",
      numerator = e1@numerator * e2,
      denominator = e1@denominator
    ),
    "/" = new(
      "ratioOfQsprays",
      numerator   = e1@numerator / e2,
      denominator = e1@denominator
    ),
    "^" = ratioOfQspraysPower(e1, e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
character_arith_ratioOfQsprays <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.ratioOfQsprays(e1) + e2,
    "-" = as.ratioOfQsprays(e1) - e2,
    "*" = new(
      "ratioOfQsprays",
      numerator = e1 * e2@numerator,
      denominator = e2@denominator
    ),
    "/" = simplifyRatioOfQsprays( # juste au cas où le num est constant
      new(
        "ratioOfQsprays",
        numerator   = e1 * e2@denominator,
        denominator = e2@numerator
      )
    ),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
gmp_arith_ratioOfQsprays <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.ratioOfQsprays(e1) + e2,
    "-" = as.ratioOfQsprays(e1) - e2,
    "*" = new(
      "ratioOfQsprays",
      numerator = e1 * e2@numerator,
      denominator = e2@denominator
    ),
    "/" = simplifyRatioOfQsprays( # juste au cas où le num est constant
      new(
        "ratioOfQsprays",
        numerator   = e1 * e2@denominator,
        denominator = e2@numerator
      )
    ),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}
numeric_arith_ratioOfQsprays <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.ratioOfQsprays(e1) + e2,
    "-" = as.ratioOfQsprays(e1) - e2,
    "*" = new(
      "ratioOfQsprays",
      numerator = e1 * e2@numerator,
      denominator = e2@denominator
    ),
    "/" = simplifyRatioOfQsprays( # juste au cas où le num est constant
      new(
        "ratioOfQsprays",
        numerator   = e1 * e2@denominator,
        denominator = e2@numerator
      )
    ),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

setMethod(
  "Arith",
  signature(e1 = "ratioOfQsprays", e2 = "ratioOfQsprays"),
  ratioOfQsprays_arith_ratioOfQsprays
)
setMethod(
  "Arith",
  signature(e1 = "ratioOfQsprays", e2 = "qspray"),
  ratioOfQsprays_arith_qspray
)
setMethod(
  "Arith",
  signature(e1 = "ratioOfQsprays", e2 = "character"),
  ratioOfQsprays_arith_character
)
setMethod(
  "Arith",
  signature(e1 = "ratioOfQsprays", e2 = "bigq"),
  ratioOfQsprays_arith_gmp
)
setMethod(
  "Arith",
  signature(e1 = "ratioOfQsprays", e2 = "bigz"),
  ratioOfQsprays_arith_gmp
)
setMethod(
  "Arith",
  signature(e1 = "qspray", e2 = "ratioOfQsprays"),
  qspray_arith_ratioOfQsprays
)
setMethod(
  "Arith",
  signature(e1 = "character", e2 = "ratioOfQsprays"),
  character_arith_ratioOfQsprays
)
setMethod(
  "Arith",
  signature(e1 = "bigq", e2 = "ratioOfQsprays"),
  gmp_arith_ratioOfQsprays
)
setMethod(
  "Arith",
  signature(e1 = "bigz", e2 = "ratioOfQsprays"),
  gmp_arith_ratioOfQsprays
)
setMethod(
  "Arith",
  signature(e1 = "ratioOfQsprays", e2 = "numeric"),
  ratioOfQsprays_arith_numeric
)
setMethod(
  "Arith",
  signature(e1 = "numeric", e2 = "ratioOfQsprays"),
  numeric_arith_ratioOfQsprays
)


setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    num <- e1@numerator*e2@denominator - e2@numerator*e1@denominator
    switch(
      .Generic,
      "==" = isQzero(num),
      "!=" = !(isQzero(num)),
      stop(gettextf(
        "Comparison operator %s not defined for ratioOfQsprays objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "qspray"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.ratioOfQsprays(e2),
      "!=" = e1 != as.ratioOfQsprays(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "character"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.ratioOfQsprays(e2),
      "!=" = e1 != as.ratioOfQsprays(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "numeric"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.ratioOfQsprays(e2),
      "!=" = e1 != as.ratioOfQsprays(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "bigz"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.ratioOfQsprays(e2),
      "!=" = e1 != as.ratioOfQsprays(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "bigq"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.ratioOfQsprays(e2),
      "!=" = e1 != as.ratioOfQsprays(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "qspray", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.ratioOfQsprays(e1) == e2,
      "!=" = as.ratioOfQsprays(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "character", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.ratioOfQsprays(e1) == e2,
      "!=" = as.ratioOfQsprays(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "numeric", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.ratioOfQsprays(e1) == e2,
      "!=" = as.ratioOfQsprays(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "bigz", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.ratioOfQsprays(e1) == e2,
      "!=" = as.ratioOfQsprays(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "bigq", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.ratioOfQsprays(e1) == e2,
      "!=" = as.ratioOfQsprays(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.",
        dQuote(.Generic)
      ))
    )
  }
)
