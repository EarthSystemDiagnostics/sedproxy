library(sedproxy)
context("ClimToProxyClim")

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


# test_that("First rep equal to single rep"){
#
#   clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
#   clim.in <- ts(clim.in, start = -39)
#
#   tpts <- round(seq(4000, 20000, length.out = 10))
#
#   set.seed(26052017)
#   PFM.slow.1 <- ClimToProxyClim(clim.signal = clim.in,
#                                 timepoints = tpts,
#                                 habitat.weights = N41.G.ruber.seasonality,
#                                 sed.acc.rate = rep(50, length(tpts)),
#                                 layer.width = 0,
#                                 sigma.meas = 0.46,
#                                 sigma.ind = 0,
#                                 n.samples = 30,
#                                 n.replicates = 1)
#   set.seed(26052017)
#   PFM.slow.2 <- ClimToProxyClim(clim.signal = clim.in,
#                                 timepoints = tpts,
#                                 habitat.weights = N41.G.ruber.seasonality,
#                                 sed.acc.rate = rep(50, length(tpts)),
#                                 layer.width = 0,
#                                 sigma.meas = 0.46,
#                                 sigma.ind = 0,
#                                 n.samples = 30,
#                                 n.replicates = 3)
#
#   expect_equal(object = dplyr::filter(data.frame(PFM.slow.1$everything), replicate == 1),
#                expected = dplyr::filter(data.frame(PFM.slow.2$everything), replicate == 1))
#
#   expect_equal(object = data.frame(PFM.slow.1$simulated.proxy),
#                expected = data.frame(PFM.slow.2$simulated.proxy))
#
# }

