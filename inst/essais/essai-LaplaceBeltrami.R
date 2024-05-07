library(ratioOfQsprays)

LaplaceBeltrami <- function(qspray, alpha) {
  n <- numberOfVariables(qspray)
  derivatives <- lapply(seq_len(n), function(i) {
    derivQspray(qspray, i)
  })
  x <- lapply(seq_len(n), qlone) # x_1, x_2, ..., x_n
  out <- 0L
  for(i in seq_len(n)) {
    out <- out + alpha * (x[[i]] * derivatives[[i]])^2
    for(j in seq_len(n)) {
      if(i != j) {
        out <- out + x[[i]]^2 * derivatives[[i]] / (x[[i]] - x[[j]])
      }
    }
  }
  # at this step, `out` is a `ratioOfQsprays` object, because of the divisions
  # by `x[[i]] - x[[j]]`; but actually its denominator is 1 because of some
  # simplifications and then we extract its numerator to get a `qspray` object
  getNumerator(out) / 2
}

( qspray <- PSFpoly(3, c(3,2,1)) )
LaplaceBeltrami(qspray) |> isPolynomial()
