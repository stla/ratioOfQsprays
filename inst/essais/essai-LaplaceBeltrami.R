library(ratioOfQsprays)

LaplaceBeltrami <- function(qspray) {
  n <- numberOfVariables(qspray)
  derivatives <- lapply(seq_len(n), function(i) {
    derivQspray(qspray, i)
  })
  x <- lapply(seq_len(n), qlone)
  out <- 0L
  for(i in seq_len(n)) {
    for(j in seq_len(n)) {
      if(i != j) {
        out <- out + x[[i]]^2 * derivatives[[i]] / (x[[i]] - x[[j]])
      }
    }
  }
  out
}

( qspray <- PSFpoly(3, c(3,2,1)) )
LaplaceBeltrami(qspray) |> isPolynomial()
