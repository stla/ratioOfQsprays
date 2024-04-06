passShowAttributes <- function(source, target) {
  lapply(c("x", "quotientBar"), function(a) {
    attr(target, a) <<- attr(source, a)
  })
  target
}

isInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && as.integer(x) == x
}

isPositiveInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isNonnegativeInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x && x != 0
}

#' @importFrom qspray numberOfVariables
#' @noRd
numberOfVariables2 <- function(roq) {
  max(
    numberOfVariables(roq@numerator),
    numberOfVariables(roq@denominator)
  )
}

#' @importFrom utils head
#' @noRd
removeTrailingZeros <- function(x) {
  n <- length(x)
  while(x[n] == 0 && n > 0L) {
    n <- n - 1L
  }
  head(x, n)
}

qspray_from_list <- function(qspray_as_list) {
  powers <- qspray_as_list[["powers"]]
  if(is.null(powers)) {
    new("qspray", powers = list(), coeffs = character(0L))
  } else {
    new("qspray", powers = powers, coeffs = qspray_as_list[["coeffs"]])
  }
}
