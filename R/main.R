getSpatialWeights <- function(pos, klist=3:9, plot=TRUE) {
  # merge many adjs
  adjList <- lapply(klist, function(k) {
    adj <- getAdj(pos, k=k)
  })
  adj <- Reduce("+", adjList) / length(adjList)
  if(plot) {
    plotNetwork(pos, adj, line.power=10)
  }
  return(adj)
}

getSpatialPatterns <- function(mat, adj, ncores=parallel::detectCores()-1) {
  results <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
    value <- mat[i,]
    ape::Moran.I(value, adj)
  }, mc.cores=ncores))
  rownames(results) <- rownames(mat)
  results <- as.data.frame(results)

  results$p.adj <- stats::p.adjust(results$p.value)
  results <- results[order(results$p.adj),]
  return(results)
}
