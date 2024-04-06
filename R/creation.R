#' @title Random 'ratioOfQsprays'
#' @description Generates a random \code{ratioOfQsprays} object.
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#' @importFrom qspray rQspray
rRatioOfQsprays <- function() {
  numerator   <- rQspray()
  denominator <- rQspray()
  while(denominator == 0L) {
    denominator <- rQspray()
  }
  numerator / denominator
}
