#' Derive spatial weights
#'
#' @description Derive spatial weights given position of cells in space
#'
#' @param pos Position matrix where each row is a cell, columns are
#'     x, y, (optionally z) coordinations
#' @param klist range of number of nearest neighbors to consider (default 3:9)
#' @param ncores Number of cores
#' @param plot Whether to plot neighbor network
#'
#' @returns A weighted adjacency matrix
#'
#' @export
#'
getSpatialWeights <- function(pos, klist=6, ncores=1, plot=FALSE, verbose=TRUE) {
  if(length(klist)>1) {
    adjList <-  BiocParallel::bplapply(seq_along(klist), function(i) {
      k <- klist[i]
      adj <- getAdj(pos, k=k)
    }, BPPARAM = BiocParallel::MulticoreParam(workers=ncores, tasks=length(klist)/ncores, progressbar=verbose))
    adj <- Reduce("+", adjList) / length(adjList)
  } else {
    adj <- getAdj(pos, k=klist)
  }
  if(plot) {
    plotNetwork(pos, adj, line.power=3)
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
#'
#' @export
#'
getSpatialPatterns <- function(mat, adj) {

  # Calculate Moran's I for each gene
  results <- getSpatialPatterns_C(as.matrix(mat), as.matrix(adj))
  colnames(results) <- c('observed', 'expected', 'sd')
  rownames(results) <- rownames(mat)
  results <- as.data.frame(results)

  # Get p-value
  # always assume we want greater autocorrelation
  pv <- 1 - pnorm(results$observed, mean = results$expected, sd = results$sd)
  results$p.value <- pv;
  # multiple testing correction
  results$p.adj <- stats::p.adjust(results$p.value)
  # order by significance
  results <- results[order(results$p.adj),]

  return(results)
}


groupSigSpatialPatterns <- function(pos, mat, d, deepSplit=0, minClusterSize=0, plot=TRUE, verbose=TRUE, ...) {

    #m <- as.matrix(mat[results.sig,])
    #m <- apply(m, 1, winsorize) # get rid of outliers
    #d <- as.dist(1-cor(m))

    hc <- hclust(d)
    if(plot) {
        par(mfrow=c(1,1))
        plot(hc)
    }
    # dynamic tree cut
    groups <- dynamicTreeCut::cutreeDynamic(dendro=hc, distM=as.matrix(d), method='hybrid', minClusterSize=minClusterSize, deepSplit=deepSplit)
    names(groups) <- hc$labels
    groups <- factor(groups)
    if(verbose) {
        message('Patterns detected:')
        print(table(groups))
    }

    if(plot) {
        par(mfrow=c(length(levels(groups)), 2), mar=rep(1,4))
    }

    prs <- lapply(levels(groups), function(g) {
        # summarize as first pc if more than 1 gene in group
        if(sum(groups==g)>1) {
            m <- winsorize(mat[results.sig[groups==g],])
            ps <- colMeans(m)
            pr <- scale(ps)[,1]
        } else {
            ps <- winsorize(mat[results.sig[groups==g],])
            pr <- scale(ps)[,1]
        }
        if(plot) {
            interpolate(pos, pr, main=paste0("Pattern ", g, " : ", sum(groups==g), " genes"), plot=TRUE, ...)
        }
        return(pr)
    })
    names(prs) <- levels(groups)

    return(list(hc=hc, groups=groups, prs=prs))
}

