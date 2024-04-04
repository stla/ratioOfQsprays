isInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && as.integer(x) == x
}

isPositiveInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isNonnegativeInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x && x != 0
}

# #' @importFrom qspray numberOfVariables
numberOfVariables <- function(qspray) {
  suppressWarnings(max(lengths(qspray@powers)))
}
numberOfVariables2 <- function(roq) {
  max(
    suppressWarnings(max(lengths(roq@numerator@powers))),
    suppressWarnings(max(lengths(roq@denominator@powers)))
    # qspray::numberOfVariables(roq@numerator),
    # qspray::numberOfVariables(roq@denominator)
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

