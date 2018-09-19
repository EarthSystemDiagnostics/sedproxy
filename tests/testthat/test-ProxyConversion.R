library(sedproxy)
context("ProxyConversion")

test_that("Forward back", {
  input <- 5:15
  expect_equal(ProxyConversion(
    proxy.value = ProxyConversion(temperature = input, calibration.type = "MgCa"),
    calibration.type = "MgCa"), input)
  expect_equal(ProxyConversion(
    proxy.value = ProxyConversion(temperature = input, calibration.type = "Uk37"),
    calibration.type = "Uk37"),  input)
})


test_that("Output shape", {
  input <- matrix(5:13, ncol = 3)
  expect_equal(dim(ProxyConversion(temperature = input,
                                   calibration.type = "MgCa")), dim(input))
  expect_equal(dim(ProxyConversion(temperature = input,
                                   calibration.type = "MgCa",
                                   point.or.sample = "sample", n = 3)), dim(input))
  expect_equal(dim(ProxyConversion(temperature = input,
                                   calibration.type = "Uk37")), dim(input))
  expect_equal(dim(ProxyConversion(temperature = input,
                                   calibration.type = "Uk37",
                                   point.or.sample = "sample", n = 3)), dim(input))
  
  
  input <- 5:15
  n.req <- 3
  expect_equal(ncol(ProxyConversion(temperature = input,
                                   calibration.type = "MgCa",
                                   point.or.sample = "sample",
                                   n = n.req)), n.req)
  
  expect_equal(ncol(ProxyConversion(temperature = input,
                                    calibration.type = "Uk37",
                                    point.or.sample = "sample",
                                    n = n.req)), n.req)
  

  input <- matrix(5:13, ncol = 3)
  
  expect_error(
    ProxyConversion(
      temperature = input,
      calibration.type = "Uk37",
      point.or.sample = "sample",
      n = 4
    ),
    "If input is matrix and point.or.sample == 'sample', n must equal ncol(input)",
    fixed = TRUE)
})

