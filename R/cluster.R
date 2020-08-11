#' Get transcriptional clusters by graph-based clustering
#'
#' @description Get transcriptional clusters by graph-based clustering
#'
#' @param pcs A matrix of principal components or gene expression to assess transcriptional similarity
#' @param k Number of nearest neighbors for clustering
#' @param method igraph method for graph-based clustering (default: cluster_louvain)
#' @param weight Whether to using weighting by transcriptional distance
#' @param verbose Verbosity
#' @param details Return detailed ouputs
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
#' com <- getClusters(pcs, k=50, verbose=TRUE)
#'
#' @export
#'
getClusters <- function (pcs, k,
                         method = igraph::cluster_louvain,
                         weight = FALSE,
                         verbose = FALSE,
                         details = FALSE) {
  if (verbose) {
    print("finding approximate nearest neighbors ...")
  }
  # nearest neighbors in PC space
  nn = RANN::nn2(pcs, k = k) ## KNN
  names(nn) <- c('idx', 'dists')

  if(weight) {
    if(verbose) {
      print('using transcriptional distance weighting')
    }
    weight <- 1/(1+ as.vector(nn$dists))
  } else {
    if(verbose) {
      print('using equal weighting')
    }
    weight <- rep(1, nrow(pcs))
  }

  if (verbose) {
    print("calculating clustering ...")
  }
  nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                     to = as.vector(nn$idx),
                     weight = weight
  )
  g <- igraph::graph_from_data_frame(nn.df, directed = FALSE)
  g <- igraph::simplify(g)
  km <- method(g)
  if (verbose) {
    mod <- igraph::modularity(km)
    if (mod < 0.3) {
      print("WARNING")
    }
    print(paste0("graph modularity: ", mod))
  }
  com <- km$membership
  names(com) <- rownames(pcs)
  com <- factor(com)
  if (verbose) {
    print("identifying cluster membership ...")
    print(table(com))
  }
  if (details) {
    return(list(com = com, mod = mod, g = g))
  }
  else {
    return(com)
  }
}

#' Get spatially informed clusters by weighting graph-based clustering with spatial information
#'
#' @description Get spatially informed clusters by weighting graph-based clustering with spatial information
#'
#' @param pcs A matrix of principal components or gene expression to assess transcriptional similarity
#' @param W Binary adjacency matrix
#' @param k Number of nearest neighbors for clustering
#' @param method igraph method for graph-based clustering (default: cluster_louvain)
#' @param verbose Verbosity
#' @param details Return detailed ouputs
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
#' W <- getSpatialNeighbors(pos, filterDist=5)
#' com <- getSpatiallyInformedClusters(pcs, W, k=50, verbose=TRUE)
#'
#' @export
#'
getSpatiallyInformedClusters <- function(pcs, W, k,
                                         alpha=1, beta=0,
                                         method = igraph::cluster_louvain,
                                         verbose=FALSE,
                                         details=FALSE) {

  if (verbose) {
    print("finding approximate nearest neighbors ...")
  }
  # nearest neighbors in PC space
  nn = RANN::nn2(pcs, k = k) ## KNN
  names(nn) <- c('idx', 'dists')

  if (verbose) {
    print("using spatial weights ...")
  }
  # create weights from binary adjacency matrix W
  nw.simple = igraph::graph_from_adjacency_matrix(W)
  nw.simple = igraph::simplify(nw.simple)
  pweight <- do.call(rbind, lapply(1:nrow(nn$idx), function(i) {
    d <- unlist(lapply(nn$idx[i,], function(j) {
      if(i==j) { return(0) }
      else {
        return(igraph::distances(nw.simple, rownames(W)[i], colnames(W)[j]))
        #return(igraph::min_cut(nw.simple, rownames(W)[i], colnames(W)[j]))
      }
    }))
    return(d)
  }))

  if (verbose) {
    print(paste0("calculating weights with alpha:", alpha, " and beta:", beta))
  }
  weight <- 1/(alpha + as.vector(pweight)) + beta

  if (verbose) {
    print("calculating clustering ...")
  }
  # graph-based clustering with spatial weights
  nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                     to = as.vector(nn$idx),
                     weight = weight
  )
  g <- igraph::graph_from_data_frame(nn.df, directed = FALSE)
  g <- igraph::simplify(g)
  km <- method(g)
  if (verbose) {
    mod <- igraph::modularity(km)
    if (mod < 0.3) {
      print("WARNING")
    }
    print(paste0("graph modularity: ", mod))
  }
  com <- km$membership
  names(com) <- rownames(pcs)
  com <- factor(com)
  if (verbose) {
    print("identifying cluster membership ...")
    print(table(com))
  }
  if (details) {
    return(list(com = com, mod = mod, g = g))
  }
  else {
    return(com)
  }
}


