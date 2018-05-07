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
#' @param permutation Whether to use a permutation-based testing. Empirical
#'     test used by default
#' @param ncores Number of cores for parallelization across genes
#' @param ... Additional parameters to send to moranPermutationTest
#'
#' @export
#'
getSpatialPatterns <- function(mat, adj, permutation=FALSE, ncores=1, verbose=TRUE, ...) {

  ## if(verbose) {
  ##   pb <- txtProgressBar(min=0, max=nrow(mat), char = ".") # use progress bar
  ## }
  results <- do.call(rbind, BiocParallel::bplapply(seq_len(nrow(mat)), function(i) {
    ## if(verbose) {
    ##   setTxtProgressBar(pb, i)
    ## }
    value <- mat[i,]
    if(permutation) {
      moranPermutationTest(value, adj, ncores=1, ...)
    } else {
      moranTest(value, adj)
    }
  }, BPPARAM = BiocParallel::MulticoreParam(workers=ncores, tasks=nrow(mat), progressbar=verbose)))
  ## if(verbose) {
  ##   close(pb)
  ## }

  rownames(results) <- rownames(mat)
  results <- as.data.frame(results)

  results$p.adj <- stats::p.adjust(results$p.value)
  results <- results[order(results$p.adj),]

  return(results)
}


groupSigSpatialPatterns <- function(pos, mat, results, alpha=0.05, k=5, plot=TRUE, verbose=TRUE, ...) {
    vi <- results$p.adj < alpha
    vi[is.na(vi)] <- FALSE
    if(verbose) {
        message(paste0('Number of significantly spatially clustered genes: ', sum(vi)))
    }
    results.sig <- rownames(results)[vi]

    m <- as.matrix(mat[results.sig,])
    m <- apply(m, 1, winsorize) # get rid of outliers
    m <- scale(m)
    d <- dist(t(m))
    hc <- hclust(d)
    # dynamic tree cut
    groups <- dynamicTreeCut::cutreeDynamic(dendro=hc, distM=as.matrix(d), method='hybrid', minClusterSize=10, deepSplit=0)
    names(groups) <- hc$labels
    groups <- factor(groups)
    if(verbose) {
        message('Patterns detected:')
        message(table(groups))
    }

    if(plot) {
        par(mfrow=c(length(levels(groups)), 2), mar=rep(1,4))
    }

    prs <- lapply(levels(groups), function(g) {
        # summarize as first pc if more than 1 gene in group
        if(sum(groups==g)>1) {
            m <- winsorize(mat[results.sig[groups==g],])
            ps <- colMeans(m)
            pc <- prcomp(m)
            pr <- pc$rotation[,1]
            # double check direction is same
            # if not, flip pc
            if(cor(ps, pr)<0) {
                pr <- -pr
            }
        } else {
            pr <- winsorize(mat[results.sig[groups==g],])
        }
        if(plot) {
            interpolate(pos, pr, main=paste0("Pattern ", g, " : ", sum(groups==g), " genes"), plot=TRUE, ...)
        }
        return(pr)
    })
    names(prs) <- levels(groups)

    return(list(groups=groups, prs=prs))
}
