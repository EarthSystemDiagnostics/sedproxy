library(sedproxy)
context("FAME weights")

test_that("FAME weights", {
  
  w.ruber <- growth_rate_l09_R("ruber", seq(0, 40) + 273.15)
  w.ruber.thresh <- growth_rate_l09_R("ruber", seq(0, 40) + 273.15, min.growth.thresh = 0.01)
  
  w.pachy_d.norm <- growth_rate_l09_R("pachy_d", seq(18, 22, 0.1) + 273.15, norm = 1)
  
  
    
  expect_type(w.ruber, "double")
  expect_length(w.ruber, 41)
  
  
  expect_type(w.ruber.thresh, "double")
  expect_length(w.ruber.thresh, 41)
  
  expect_true(all(w.ruber.thresh[c(1:13, 34:41)] == 0))
  
  expect_gt(max(w.pachy_d.norm), 0.9)
  
})


