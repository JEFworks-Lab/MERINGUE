#' Derive spatial weights
#'
#' @description Derive spatial weights given position of cells in space
#'
#' @param pos Position matrix where each row is a cell, columns are
#'     x, y, (optionally z) coordinations
#' @param klist range of number of nearest neighbors to consider (default 3:9)
#' @param plot Whether to plot neighbor network
#'
#' @returns A weighted adjacency matrix
#'
#' @export
#'
getSpatialWeights <- function(pos, klist=3:9, plot=FALSE) {
  adjList <- lapply(klist, function(k) {
    adj <- getAdj(pos, k=k)
  })
  adj <- Reduce("+", adjList) / length(adjList)
  if(plot) {
    plotNetwork(pos, adj, line.power=10)
  }
  return(adj)
}


#' Identify spatial clusters
#'
#' @description Identify spatially clustered genes using Moran's I
#'
#' @param mat Gene expression matrix. Must be normalized such that correlations
#'     will not be driven by technical artifacts.
#' @param adj Spatial weights such as a weighted adjacency matrix
#' @param permutation Whether to use a permutation-based testing. Empirical
#'     test used by default
#' @param ncores Number of cores for parallelization across genes
#' @param ... Additional parameters to send to moranPermutationTest
#'
#' @export
#'
getSpatialPatterns <- function(mat, adj, permutation=FALSE, ncores=parallel::detectCores()-1, ...) {
  results <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
    value <- mat[i,]
    if(permutation) {
      moranPermutationTest(value, adj, ncores=1, ...)
    } else {
      moranTest(value, adj)
    }
  }, mc.cores=ncores))
  rownames(results) <- rownames(mat)
  results <- as.data.frame(results)

  results$p.adj <- stats::p.adjust(results$p.value)
  results <- results[order(results$p.adj),]
  return(results)
}


groupSpatialPatterns <- function(pos, mat) {

}



