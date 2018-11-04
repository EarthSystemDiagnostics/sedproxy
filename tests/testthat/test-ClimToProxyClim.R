library(sedproxy)
context("ClimToProxyClim")

test_that("Fast = Slow", {

  clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
  clim.in <- ts(clim.in, start = -39)

  tpts <- round(seq(40, 22000, length.out = 100))

  set.seed(26052017)
  PFM.slow <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              habitat.weights = N41.G.ruber.seasonality,
                              sed.acc.rate = rep(50, length(tpts)),
                              layer.width = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)
  set.seed(26052017)
  PFM.fast <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              habitat.weights = N41.G.ruber.seasonality,
                              sed.acc.rate = 50,
                              layer.width = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)

  set.seed(26052017)
  PFM.slow.inf <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              habitat.weights = N41.G.ruber.seasonality,
                              sed.acc.rate = rep(50, length(tpts)),
                              layer.width = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = Inf,
                              n.replicates = 1)
  set.seed(26052017)
  PFM.fast.inf <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              habitat.weights = N41.G.ruber.seasonality,
                              sed.acc.rate = 50,
                              layer.width = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = Inf,
                              n.replicates = 1)

  expect_equal(object = data.frame(PFM.fast$simulated.proxy),
                         expected = data.frame(PFM.slow$simulated.proxy))
  expect_equal(object = data.frame(PFM.fast$everything),
               expected = data.frame(PFM.slow$everything))

  expect_equal(object = data.frame(PFM.fast.inf$simulated.proxy),
               expected = data.frame(PFM.slow.inf$simulated.proxy))
  expect_equal(object = data.frame(PFM.fast.inf$everything),
               expected = data.frame(PFM.slow.inf$everything))
})


test_that("Slow = Cached slow", {
  load("data/PFM.cache.Rdata")

  clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
  clim.in <- ts(clim.in, start = -39)

  tpts <- round(seq(4000, 20000, length.out = 10))

  set.seed(26052017)
  PFM.slow <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              habitat.weights = N41.G.ruber.seasonality,
                              sed.acc.rate = rep(50, length(tpts)),
                              layer.width = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)
  #PFM.cache <- PFM.slow
  #save(PFM.cache, file = "tests/testthat/data/PFM.cache.Rdata")

  expect_equal(object = data.frame(PFM.cache$simulated.proxy),
               expected = data.frame(PFM.slow$simulated.proxy))

})

