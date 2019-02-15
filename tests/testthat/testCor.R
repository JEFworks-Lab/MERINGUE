library(testthat)

test_that(context("Simulate inter cell-type spatial cross-correlation tests"), {
  library(MERingue)

  # Simulate data
  set.seed(0)
  N <- 100
  pos <- cbind(rnorm(N), rnorm(N))
  rownames(pos) <- paste0('cell', 1:N)
  colnames(pos) <- c('x', 'y')
  ctA <- sample(rownames(pos), N/2)
  ctB <- setdiff(rownames(pos), ctA)
  gexpA <- pos[,2]
  gexpA[ctB] <- 0
  gexpB <- pos[,2]
  gexpB[ctA] <- 0
  #plotEmbedding(pos, col=gexpA)
  #plotEmbedding(pos, col=gexpB)

  weight <- getMnn(ctA, ctB, pos, k=6)
  plotNetwork(pos, weight)
  points(pos[ctA,], col='orange', pch=16)
  points(pos[ctB,], col='green', pch=16)

  cor <- cor(gexpA, gexpB)
  scor <- interCellTypeSpatialCrossCor(gexpA, gexpB, ctA, ctB, weight)

  # Test plotting
  plotEmbedding(pos, col=gexpA)
  plotEmbedding(pos, col=gexpB)
  invisible(interpolate(pos, gexpA))
  invisible(interpolate(pos, gexpB))
  plotInterCellTypeSpatialCrossCor(gexpA, gexpB, ctA, ctB, weight)

  # toroidal shift model
  set.seed(0)
  posr <- MERingue:::rtorShift(pos)
  plotEmbedding(posr, col=gexpA)
  plotEmbedding(posr, col=gexpB)
  weightr <- getMnn(ctA, ctB, posr, k=6)
  plotNetwork(posr, weightr)
  points(posr[ctA,], col='orange', pch=16)
  points(posr[ctB,], col='green', pch=16)
  scorr <- interCellTypeSpatialCrossCor(gexpA, gexpB, ctA, ctB, weightr)

  # random label model
  set.seed(0)
  gexpAr <- gexpA
  gexpBr <- gexpB
  names(gexpAr) <- sample(names(gexpA), length(gexpA), replace = FALSE)
  names(gexpBr) <- sample(names(gexpB), length(gexpB), replace = FALSE)
  plotEmbedding(pos, col=gexpAr)
  plotEmbedding(pos, col=gexpBr)
  scorr2 <- interCellTypeSpatialCrossCor(gexpAr, gexpBr, ctA, ctB, weight)

  expect_equal(scor > cor, TRUE) # regular correlation
  expect_equal(scor > scorr, TRUE) # random toroidal shift
  expect_equal(scor > scorr2, TRUE) # random label
})
