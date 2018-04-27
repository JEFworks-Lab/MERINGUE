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


groupSigSpatialPatterns <- function(pos, mat, results, alpha=0.05, k=5, verbose=TRUE) {
  vi <- results$p.adj < alpha
  if(verbose) {
    print(table(vi))
  }
  results.sig <- rownames(results)[vi]

  cv <- cor(t(as.matrix(mat[results.sig,])))
  hc <- hclust(as.dist(1-cv))
  heatmap(cv[hc$labels, hc$labels], Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorRampPalette(c('blue', 'white', 'red'))(100), labRow=NA, labCol=NA)
  groups <- cutree(hc, k)
  groups <- factor(groups)

  par(mfrow=c(length(levels(groups)), 2), mar=rep(1,4))
  prs <- lapply(levels(groups), function(g) {
    pc <- prcomp(mat[results.sig[groups==g],])
    pr <- pc$rotation[,1]
    interpolate(pos, pr)
    return(pr)
  })
  names(prs) <- levels(groups)

  return(list(groups=groups, prs=prs))
}



