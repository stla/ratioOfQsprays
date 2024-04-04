test_that("simple expansion", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  expect_true((roq1 - roq2)^2 == roq1^2 - 2*roq1*roq2 + roq2^2)
})

test_that("commutativity", {
  roq1 <- ROQ1()
  roq2 <- ROQ2()
  expect_true(roq1*roq2 == roq2*roq1)  
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
  expect_true(roq1*(roq2+roq3) == roq1*roq2 + roq1*roq3)
})