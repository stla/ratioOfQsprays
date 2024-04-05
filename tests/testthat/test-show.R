test_that("simple expansion", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  expect_true(Print((roq1 - roq2)^2) == Print(roq1^2 - 2*roq1*roq2 + roq2^2))
})

test_that("commutativity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  expect_true(Print(roq1*roq2*3) == Print(3*roq2*roq1))
  x <- qlone(1)
  expect_true(Print(roq1*roq2*x) == Print(x*roq2*roq1))
})

test_that("associativity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(Print(roq1*(roq2*roq3)) == Print((roq1*roq2)*roq3))
})

test_that("distributivity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(
    Print(roq1*(roq2+roq3+3)) == Print(roq1*roq2 + roq1*roq3 + roq1*3)
  )
})

test_that("equality with scalar", {
  roq1 <- ROQ1()
  expect_true(Print(4*roq1/(2*roq1)) == "[2*x^()] ")
  expect_true(Print(2*roq1/(4*roq1)) == "[1/2*x^()] ")
})

test_that("division", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(Print((roq1/roq2) * roq3) == Print((roq1*roq3) / roq2))
})

test_that("power", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(Print((roq1/roq2*roq3)^3) == Print(roq1^3/roq2^3*roq3^3))
})

test_that("arithmetic between qsprays and ratioOfQsprays", {
  x <- qlone(1)
  y <- qlone(2)
  z <- qlone(3)
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(
    Print("3/2"+(x*roq1)/(y*roq2)+z) == Print(z+((x/y)*(roq1/roq2))+"3/2")
  )
})

test_that("equality between qspray and ratioOfQsprays", {
  x <- qlone(1)
  y <- qlone(2)
  z <- qlone(3)
  expect_true(Print((x^2-y^2)/(x+y)) == "[x^(1) - x^(0, 1)] ")
})

test_that("equality between scalar and ratioOfQsprays", {
  roq1 <- ROQ1()
  expect_true(Print((3*roq1)/(roq1*6)) == "[1/2*x^()] ")
})
