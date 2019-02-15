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

  cor <- cor(gexpA, gexpB)
  scor <- interCellTypeSpatialCrossCor(gexpA, gexpB, ctA, ctB, weight)

  expect_equal(scor > cor, TRUE)
})
