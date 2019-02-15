library(testthat)

test_that(context("Test winsorization"), {

  library(MERingue)
  x <- rnorm(100,0,1)
  x <- c(x, 10)
  xw <- winsorize(x, 0.01)
  expect_equal(max(x[1:100]), xw[101])
})

test_that(context("Test differential expression"), {

  library(MERingue)

  set.seed(0)
  G <- 2
  N <- 30
  M <- 1000
  initmean <- 5
  initvar <- 10
  mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)
  mat <- abs(mat)
  rownames(mat) <- paste0('gene', 1:M)
  colnames(mat) <- paste0('cell', 1:(N*G))
  group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
  names(group) <- colnames(mat)
  #heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

  set.seed(0)
  upreg <- 100
  upregvar <- 10
  ng <- 100

  diff <- lapply(1:G, function(x) {
    diff <- rownames(mat)[(((x-1)*ng)+1):(((x-1)*ng)+ng)]
    mat[diff, group==paste0('group', x)] <<- mat[diff, group==paste0('group', x)] + rnorm(ng, upreg, upregvar)
    return(diff)
  })
  names(diff) <- paste0('group', 1:G)

  mat <- round(mat)
  #heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

  dg <- getDifferentialGenes(mat, group)
  dg.sig <- lapply(dg, function(x) {
    na.omit(rownames(x)[which(x$Z>3)])
  })

  expect_equal(length(intersect(dg.sig[[1]], diff[[1]])), upreg)
  expect_equal(length(intersect(dg.sig[[2]], diff[[2]])), upreg)
})
