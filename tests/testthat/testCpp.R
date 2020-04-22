library(testthat)

test_that(context("Moran's I C++ functions compiled and works as expected"), {
  library(MERingue)
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  mat <- normalizeCounts(cd, verbose=FALSE)
  w <- getSpatialNeighbors(pos)

  set.seed(0)
  is <- sample(1:nrow(mat), 100)

  start_time <- Sys.time()
  moran <- do.call(cbind, lapply(is, function(i) {
    moranTest_DEPRECATED(mat[i,], w)
  }))
  end_time <- Sys.time()
  moranTime <- end_time - start_time

  start_time <- Sys.time()
  moranC <- do.call(cbind, lapply(is, function(i) {
    moranTest(mat[i,], w)
  }))
  end_time <- Sys.time()
  moranCTime <- end_time - start_time

  #expect_equal(moranCTime < moranTime, TRUE)
  expect_equal(all.equal(moran[1,], moranC[1,]), TRUE)
  expect_equal(all.equal(moran[2,], moranC[2,]), TRUE)
  expect_equal(all.equal(moran[3,], moranC[3,]), TRUE)
  expect_equal(all.equal(moran[4,], moranC[4,]), TRUE)

  #x <- moranP <- moranPermutationTest(mat[is[1],], w)
  #y <- moranTest(mat[is[1],], w)
  #all.equal(x[1], y[1])
  #all.equal(x[4], y[4])
})

test_that(context("Spatial cross correlation C++ functions compiled and works as expected"), {
  library(MERingue)
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  mat <- normalizeCounts(cd, verbose=FALSE)
  w <- getSpatialNeighbors(pos)

  set.seed(1)
  is <- sample(1:nrow(mat), 100)

  scc <- spatialCrossCorMatrix(as.matrix(mat[is,]), w)
  moran <- sapply(is, function(i) {
    moranTest(mat[i,], w)[1]
  })
  names(moran) <- rownames(mat)[is]

  expect_equal(all.equal(moran, diag(scc)), TRUE)
  expect_equal(all.equal(spatialCrossCor(mat[is[1],], mat[is[2],], w), scc[1,2]), TRUE)
  expect_equal(all.equal(spatialCrossCor(mat[is[2],], mat[is[3],], w), scc[2,3]), TRUE)
  expect_equal(all.equal(spatialCrossCor(mat[is[20],], mat[is[10],], w), scc[20,10]), TRUE)
})

test_that(context("LISA works as expected"), {
  library(MERingue)
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  mat <- normalizeCounts(cd, verbose=FALSE)
  w <- getSpatialNeighbors(pos)

  set.seed(0)
  is <- sample(1:nrow(mat), 10)

  moranC <- sapply(is, function(i) {
    moranTest(mat[i,], w)[1]
  })
  mLisa <- sapply(is, function(i) {
    mean(lisaTest(mat[i,], w)[,1])
  })
  names(moranC) <- names(mLisa) <- rownames(mat)[is]

  expect_equal(all.equal(moranC, mLisa), TRUE)
})

test_that(context("getSpatialPatterns works"), {
  library(MERingue)
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  mat <- normalizeCounts(cd, log=FALSE, verbose=FALSE)
  w <- getSpatialNeighbors(pos)

  set.seed(0)
  is <- rownames(mat)[sample(1:nrow(mat), 10)]

  # gold standard
  start_time <- Sys.time()
  I1 <- do.call(rbind, lapply(is, function(g) { moranTest(mat[g,], w) }))
  rownames(I1) <- is
  end_time <- Sys.time()
  moranTime <- end_time - start_time

  # test
  start_time <- Sys.time()
  I2 <- getSpatialPatterns(mat[is,], w)
  end_time <- Sys.time()
  moranCTime <- end_time - start_time

  expect_equal(moranCTime < moranTime, TRUE)
  expect_equal(all.equal(as.numeric(I1[,1]), as.numeric(I2[,1])), TRUE)

})
