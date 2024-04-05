test_that("simple expansion", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  expect_true((roq1 - roq2)^2 == roq1^2 - 2*roq1*roq2 + roq2^2)
})

test_that("commutativity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  expect_true(roq1*roq2*3 == 3*roq2*roq1)
  x <- qlone(1)
  expect_true(roq1*roq2*x == x*roq2*roq1)
})

test_that("associativity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(roq1*(roq2*roq3) == (roq1*roq2)*roq3)
})

test_that("distributivity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true(roq1*(roq2+roq3+3) == roq1*roq2 + roq1*roq3 + roq1*3)
})

test_that("equality", {
  roq1 <- ROQ1()
  expect_true(4*roq1/(2*roq1) == 2L)
  expect_true(4*roq1/(2*roq1) == as.character(2L))
  expect_true(4*roq1/(2*roq1) == gmp::as.bigq(2L))
  expect_true(4*roq1/(2*roq1) == as.qspray(2L))
  expect_true(4*roq1/(2*roq1) == as.ratioOfQsprays(2L))
})

test_that("division", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  roq3 <- ROQ3()
  expect_true((roq1/roq2) * roq3 == (roq1*roq3) / roq2)
})
