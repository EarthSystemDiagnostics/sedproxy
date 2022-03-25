library(sedproxy)
context("ClimToProxyClim")

# Fast = Slow if no mixed layer points ------
test_that("Fast = Slow if no mixed layer points", {

  clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
  clim.in <- ts(clim.in, start = -39)

  tpts <- round(seq(400, 22000, length.out = 100))

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

# zero mixing -----
test_that("zero mixing works", {

  clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
  #clim.in <- cbind(1:22040)
  clim.in <- ts(clim.in, start = 1)

  tpts <- round(seq(4000, 20000, length.out = 10))

  set.seed(26052017)
  PFM.slow <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              #habitat.weights = rep(1/12, 12),
                              sed.acc.rate = rep(50, length(tpts)),
                              layer.width = 0,
                              bio.depth = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)

  PFM.fast <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              #habitat.weights = rep(1/12, 12),
                              sed.acc.rate = 50,
                              layer.width = 0,
                              bio.depth = 0,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)

  expect_equal(object = PFM.slow$simulated.proxy$proxy.bt,
               expected = rowMeans(clim.in)[tpts])

  expect_equal(object = PFM.fast$simulated.proxy$proxy.bt,
               expected = rowMeans(clim.in)[tpts])
})

# negative times --------
test_that("negative times work -fast", {

  set.seed(26052017)
  clim.in <- matrix(rnorm(101*12), ncol = 12)
  clim.in <- ts(clim.in, start = -49)

  tpts <- c(1, 10, 100, 200)-51

  PFM <- ClimToProxyClim(clim.signal = clim.in,
                         timepoints = tpts,
                         calibration.type = "identity",
                         sed.acc.rate = 50,
                         bio.depth = 0, layer.width = 0,
                         n.samples = 30, n.replicates = 10, n.bd = 3,
                         top.of.core = -49)

  PFM$simulated.proxy$clim.signal.ann
  rowMeans(clim.in[time(clim.in)%in%tpts,])


  expect_equal(object = PFM$simulated.proxy$clim.signal.ann,
               expected = rowMeans(clim.in[time(clim.in)%in%tpts,]))

})

test_that("negative times work - slow", {
  
  set.seed(26052017)
  clim.in <- matrix(rnorm(101*12), ncol = 12)
  clim.in <- ts(clim.in, start = -49)
  
  tpts <- c(1, 10, 100, 200)-51
  
  PFM <- ClimToProxyClim(clim.signal = clim.in,
                         timepoints = tpts,
                         calibration.type = "identity",
                         sed.acc.rate = rep(50, length(tpts)),
                         bio.depth = 0, layer.width = 0,
                         n.samples = 30, n.replicates = 10, n.bd = 3,
                         top.of.core = -49)
  
  PFM$simulated.proxy$clim.signal.ann
  rowMeans(clim.in[time(clim.in)%in%tpts,])
  
  
  expect_equal(object = PFM$simulated.proxy$clim.signal.ann,
               expected = rowMeans(clim.in[time(clim.in)%in%tpts,]))
  
})

# mean bias repwise -----
test_that("meas.bias applied repwise", {
  
  clim.in <- ts(matrix(rep(1, 12*1e04), ncol = 12))
  PFM.fast <- ClimToProxyClim(clim.in,
                         timepoints = seq(1000, 6000, 500),
                         meas.bias = 3,
                         n.replicates = 3,
                        )
  
  #PlotPFMs(PFM.fast, max.replicates = 30)
  
  PFM.sub.fast <- PFM.fast$everything[PFM.fast$everything$stage == "simulated.proxy", ]
  
  PFM.slow <- ClimToProxyClim(clim.in,
                         timepoints = seq(1000, 6000, 500),
                         sed.acc.rate = rep(50, 11),
                         meas.bias = 3,
                         n.replicates = 3,
                         n.samples = 30)
  #PlotPFMs(PFM.slow)
  
  PFM.sub.slow <- PFM.slow$everything[PFM.slow$everything$stage == "simulated.proxy", ]
  
  
  PFM.mult <- ClimToProxyClim(clim.in,
                              timepoints = seq(1000, 6000, 500),
                              sed.acc.rate = rep(50, 11),
                              meas.bias = log(3),
                              n.replicates = 3,
                              n.samples = 30,
                              noise.type = "multiplicative")
  #PlotPFMs(PFM.mult)
  
  PFM.sub.mult <- PFM.mult$everything[PFM.mult$everything$stage == "simulated.proxy", ]
  
  
  expect_equal(as.numeric(tapply(PFM.sub.fast$value, PFM.sub.fast$replicate, sd)),
               c(0,0,0))
  
  expect_equal(as.numeric(tapply(PFM.sub.slow$value, PFM.sub.slow$replicate, sd)),
               c(0,0,0))
  
  expect_equal(as.numeric(tapply(PFM.sub.mult$value, PFM.sub.mult$replicate, sd)),
               c(0,0,0))
  
  })


# mixed layer -----

test_that("mixed layer modelled", {
  
  clim.in <- cbind(1:22040)
  clim.in <- ts(clim.in, start = 1)
  
  tpts <- round(seq(10, 5000, length.out = 10))
  
  set.seed(26052017)
  PFM.slow <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              #habitat.weights = rep(1/12, 12),
                              sed.acc.rate = rep(2, length(tpts)),
                              layer.width = 0,
                              bio.depth = 10,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)
  set.seed(26052017)
  PFM.fast <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              #habitat.weights = rep(1/12, 12),
                              sed.acc.rate = 2,
                              layer.width = 0,
                              bio.depth = 10,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)
  
  
  expect_equal(object = max(diff(PFM.slow$simulated.proxy$proxy.bt)),
               expected = 0)
  
  expect_equal(object = max(diff(PFM.fast$simulated.proxy$proxy.bt)),
               expected = 0)
})

test_that("mixed layer, varying sed.acc.rate", {
  
  ## Does not work with varying S
  
  clim.in <- cbind(1:22040)
  clim.in <- ts(clim.in, start = 1)
  
  tpts <- round(seq(10, 5000, length.out = 10))
  s.rates <- c(6,6,6,2,2,2, 2,2,2,2)
  #s.rates <- rep(2.2, 10)
  
  set.seed(26052017)
  PFM.slow <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              #habitat.weights = rep(1/12, 12),
                              sed.acc.rate = s.rates,
                              layer.width = 0,
                              bio.depth = 10,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)
  
  # PlotPFMs(PFM.slow) +
  #  ggplot2::geom_vline(xintercept = 1000* 10 / min(s.rates[1:4]))
  # 
  expect_equal(object = max(diff(PFM.slow$simulated.proxy$proxy.bt)),
               expected = 0)
})

# example from paper --------
test_that("example from paper works", {
  
  load("data/PFM.ex.Rdata")
  PFM.cache <- PFM.ex
  clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
  clim.in <- ts(clim.in, start = 1)
  
  
  tpts <- N41.proxy$Published.age
  s.rates <- N41.proxy$Sed.acc.rate.cm.ka
  
  
  set.seed(26052017)
  PFM.ex <- ClimToProxyClim(clim.signal = clim.in,
                              timepoints = tpts,
                              habitat.weights = N41.G.ruber.seasonality,
                              sed.acc.rate = s.rates,
                              layer.width = 1,
                              bio.depth = 10,
                              sigma.meas = 0.46,
                              sigma.ind = 0,
                              n.samples = 30,
                              n.replicates = 1)
  
  
  #save(PFM.ex, file = "tests/testthat/data/PFM.ex.Rdata")
  
  expect_equal(object = data.frame(PFM.cache$simulated.proxy),
               expected = data.frame(PFM.ex$simulated.proxy))
  
  p <- PlotPFMs(PFM.ex) 
  testthat::expect_equal(class(p), c("gg", "ggplot"))
  
  
})



