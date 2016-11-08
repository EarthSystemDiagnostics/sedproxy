library(sedproxy)
context("ImpulseResponse")

test_that("pz = 0 for depths >= d", {
  expect_equal(ImpulseResponse(z = 10.01, d = 10), 0)
})
