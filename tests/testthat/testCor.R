library(testthat)

test_that(context("Simulate inter cell-type spatial cross-correlation tests"), {
  library(MERINGUE)

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

  #weight <- getSpatialNeighbors(ctA, ctB, pos, k=6)
  weight <- getSpatialNeighbors(pos, filterDist = 1)
  weightIc <- getInterCellTypeWeight(ctA, ctB,
                                     weight, pos,
                                     plot=TRUE,
                                     main='Adjacency Weight Matrix\nBetween Cell-Types')

  cor <- cor(gexpA, gexpB)
  scor <- spatialCrossCor(gexpA, gexpB, weightIc)

  # Test plotting
  plotEmbedding(pos, col=gexpA)
  plotEmbedding(pos, col=gexpB)
  invisible(interpolate(pos, gexpA))
  invisible(interpolate(pos, gexpB))
  plotInterCellTypeSpatialCrossCor(gexpA, gexpB, ctA, ctB, weight)

  # random label model
  nullmodel <- sapply(1:1000, function(i) {
    set.seed(i)
    gexpAr <- gexpA
    gexpBr <- gexpB
    names(gexpAr) <- sample(names(gexpA), length(gexpA), replace = FALSE)
    names(gexpBr) <- sample(names(gexpB), length(gexpB), replace = FALSE)
    #plotEmbedding(pos, col=gexpAr)
    #plotEmbedding(pos, col=gexpBr)
    scorr <- spatialCrossCor(gexpAr, gexpBr, weightIc)
    #scorr <- interCellTypeSpatialCrossCor(gexpAr, gexpBr, ctA, ctB, weight)
    return(scorr)
  })
  hist(nullmodel)
  abline(v = scor, col='red')
  scorr <- mean(nullmodel)

  expect_equal(scor > cor, TRUE) # regular correlation
  expect_equal(scor > scorr, TRUE) # random label
})
