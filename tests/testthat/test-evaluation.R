test_that("evaluation", {
  library(gmp)
  f <- function(x, y, z) {
    (x + 2*y^2 - 3*z^3) / (3*x^3 - 2*y^2 + z)
  }
  roq <- f(qlone(1), qlone(2), qlone(3))
  x <- as.bigq("2")
  y <- as.bigq("3/2")
  z <- as.bigq("4/3")
  expect_true(
    evalRatioOfQsprays(roq, c(x, y, z)) == f(x, y, z)
  )
})
