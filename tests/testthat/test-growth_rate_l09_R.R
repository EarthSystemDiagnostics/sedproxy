library(sedproxy)
context("Growth rate function")

test_that("Growth rate function", {
  
  w.ruber <- ForamGrowthfT("ruber", seq(0, 40) + 273.15)
  w.ruber.thresh <- ForamGrowthfT("ruber", seq(0, 40) + 273.15,
                                  min.growth.thresh = 0.01)
  
  w.pachy_d.norm <- ForamGrowthfT("pachy_d", seq(18, 22, 0.1) + 273.15,
                                  norm = 1)
  
  
    
  expect_type(w.ruber, "double")
  expect_length(w.ruber, 41)
  
  
  expect_type(w.ruber.thresh, "double")
  expect_length(w.ruber.thresh, 41)
  
  expect_true(all(w.ruber.thresh[c(1:13, 34:41)] == 0))
  
  expect_gt(max(w.pachy_d.norm), 0.9)
  
  
  
  expected1.cache <- readRDS("data/expected1.Rdata")
  expected.norm.cache <- readRDS("data/expected.norm.Rdata")
  
  taxa <- c("sacculifer", "bulloides", "pachy_d", "siphonifera", "universa", 
          "pachy_s", "dutertrei", "ruber")

  expected1 <- lapply(taxa, function(x)
    ForamGrowthfT(x, (-5:35) + 273.15, norm = FALSE, min.growth.thresh = 0))
  
  expected.norm <- lapply(taxa, function(x)
    ForamGrowthfT(x, (-5:35) + 273.15, norm = TRUE, min.growth.thresh = 0))

  #saveRDS(expected1, "data/expected1.Rdata")
  #saveRDS(expected.norm, "data/expected.norm.Rdata")
  
  expect_equal(expected1, expected1.cache)
  expect_equal(expected.norm, expected.norm.cache)
  
  # matrix in matrix out
  
  TM <- matrix(11:16, ncol = 2)
  
  GM <- ForamGrowthfT("sacculifer", TM + 273.15)
  
  expect_equal(ncol(GM), 2)
  expect_equal(nrow(GM), 3)
  expect_true(is.matrix(GM))
  
})


