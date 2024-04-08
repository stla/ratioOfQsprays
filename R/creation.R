#' @title Random 'ratioOfQsprays'
#' @description Generates a random \code{ratioOfQsprays} object.
#'
#' @param allow.zero Boolean, whether to allow to get a null 
#'   \code{ratioOfQsprays}
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#' @importFrom qspray rQspray
rRatioOfQsprays <- function(allow.zero = TRUE) {
  numerator   <- rQspray()
  if(!allow.zero) {
    while(numerator ==0L) {
      numerator <- rQspray()
    }
  }
  denominator <- rQspray()
  while(denominator == 0L) {
    denominator <- rQspray()
  }
  numerator / denominator
}
