## random gene expression
getMinPercentCells <- function(weight, mat, alpha=0.05, M=10000, seed=0, plot=TRUE) {

  ## show shuffle
  set.seed(seed)
  rand <- mat[,sample(ncol(mat))]
  colnames(rand) <- colnames(mat)

  ## assess significance of randomly permuted genes
  I <- getSpatialPatterns(rand, weight)
  falsePositives <- rownames(I)[I$p.adj < alpha]
  print(falsePositives)

  ## get lisa
  lisa <- sapply(falsePositives, function(g) {
    gexp <- rand[g,]
    Ii <- lisaTest(gexp, weight)
    lisa <- Ii$p.value
    names(lisa) <- rownames(Ii)
    sum(lisa < alpha)/length(lisa) ## percent significant
  })

  ## determine threshold to remove false positives
  mpc <- quantile(lisa, 1-alpha)
  #mpc
  #print(sum(lisa > mpc)/length(lisa))

  if(plot) {
    par(mfrow=c(1,1), mar=rep(5,4))
    plot(lisa, -log10(I[falsePositives,]$p.adj))
    abline(h=-log10(alpha), col='red')
    abline(v=mpc, col='red')
  }

  return(mpc)
}
