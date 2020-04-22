#' Get spatially informed clusters by weighting graph-based clustering with spatial information
#'
#' @description Get spatially informed clusters by weighting graph-based clustering with spatial information
#'
#' @param pcs A matrix of principal components or gene expression to assess transcriptional similarity
#' @param W Binary adjacency matrix
#' @param k Number of nearest neighbors for clustering
#'
#' @return Factor of cluster annotations
#'
#' @examples
#' # simulate 3 spatially but 2 transcriptionally distinct groups
#' N <- 300
#' M <- 30
#' # Three spatially distinct groups
#' pos1 <- cbind(rnorm(N/3), rnorm(N/3))
#' pos2 <- cbind(rnorm(N/3, 10), rnorm(N/3))
#' pos3 <- cbind(rnorm(N/3, 10), rnorm(N/3, 10))
#' pos <- rbind(rbind(pos1, pos2), pos3)
#' group <- c(rep(1, N/3), rep(2, N/3), rep(3, N/3))
#' names(group) <- rownames(pos) <- paste0('cell', 1:N)
#' # But two are transcriptionally identical
#' pcs12 <- matrix(rnorm(N*2/3*M), N*2/3, M)
#' pcs3 <- matrix(rnorm(N*1/3*M, 10), N*1/3, M)
#' pcs <- rbind(pcs12, pcs3)
#' pcs <- cbind(pcs, abs(10-pcs))
#' colnames(pcs) <- paste0('PC:', 1:ncol(pcs))
#' rownames(pcs) <- rownames(pos)
#' com <- getSpatiallyInformedClusters(pcs, W)
#'
#' @export
#'
getSpatiallyInformedClusters <- function(pcs, W, k=30) {
  # nearest neighbros in PC space
  nn = RANN::nn2(as.matrix(pcs), k = k) ## KNN
  names(nn) <- c('idx', 'dists')

  # create weights from binary adjacency matrix W
  nw.simple = igraph::graph_from_adjacency_matrix(W)
  nw.simple = igraph::simplify(nw.simple)
  pweight <- do.call(rbind, lapply(1:nrow(nn$idx), function(i) {
    #d <- as.matrix(dist(pos[nn$idx[i,],]))
    #d[1,]
    d <- unlist(lapply(nn$idx[i,], function(j) {
      if(i==j) { return(0) }
      else {
        return(igraph::distances(nw.simple, rownames(W)[i], colnames(W)[j]))
        #return(igraph::min_cut(nw.simple, rownames(W)[i], colnames(W)[j]))
      }
    }))
    return(d)
  }))

  # graph-based clustering with spatial weights
  nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                     to = as.vector(nn$idx),
                     weight = 1/(1 + as.vector(pweight))
                     #weight = as.vector(pweight)
  )
  nw.norm = igraph::graph_from_data_frame(nn.df, directed = FALSE)
  nw.norm = igraph::simplify(nw.norm)
  lc.norm = igraph::cluster_louvain(nw.norm)

  com = as.factor(igraph::membership(lc.norm))
  names(com) <- rownames(pcs)
  return(com)
}

