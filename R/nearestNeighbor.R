#' K nearest neighbors
#'
#' @description K nearest neighbors in space based on position
#'
#' @param pos Position
#' @param k Number of nearest neighbors
#'
#' @return Boolean matrix where value = 1 if two cells are considered adjacency ie. neighbors, else 0
#'
#' @export
#'
getKnn <- function(pos, k) {
  ## nearest neighbors include self so add 1
  knn <- RANN::nn2(pos, k=k+1)[[1]]
  knn <- knn[, -1]
  ## convert to adjacency matrix
  adj <- matrix(0, nrow(pos), nrow(pos))
  rownames(adj) <- colnames(adj) <- rownames(pos)
  invisible(lapply(seq_len(nrow(pos)), function(i) {
    adj[i,rownames(pos)[knn[i,]]] <<- 1
  }))
  return(adj)
}


#' Mutual nearest neighbors
#'
#' @description Mutual nearest neighbors in space between group A and B based on position
#'
#' @param ctA vector of cell names in group A
#' @param ctB vector of cell names in group B
#' @param pos Position
#' @param k Number of mutual nearest neighbors
#'
#' @return Boolean matrix where value = 1 if two cells are considered adjacency ie. neighbors, else 0
#'
#' @export
#'
getMnn <- function(ctA, ctB, pos, k) {
    # ctB's nearest ctA neighbors
    knnB <- RANN::nn2(pos[ctA,], pos[ctB,], k=k)[[1]]
    knnB.named <- matrix(ctA[knnB], nrow=nrow(knnB), ncol=ncol(knnB))
    rownames(knnB.named) <- ctB
    # ctA's nearest ctB neighbors
    knnA <- RANN::nn2(pos[ctB,], pos[ctA,], k=k)[[1]]
    knnA.named <- matrix(ctB[knnA], nrow=nrow(knnA), ncol=ncol(knnA))
    rownames(knnA.named) <- ctA

    # mutual nearest neighbors
    ## mnnB <- lapply(1:nrow(knnB.named), function(i) {
    ##     vi <- sapply(1:k, function(j) {
    ##         rownames(knnB.named)[i] %in% knnA.named[knnB.named[i,j],]
    ##     })
    ##     knnB.named[i,][vi]
    ## })
    ## names(mnnB) <- rownames(knnB.named)
    mnnA <- lapply(seq_len(nrow(knnA.named)), function(i) {
        vi <- sapply(seq_len(k), function(j) {
            rownames(knnA.named)[i] %in% knnB.named[knnA.named[i,j],]
        })
        knnA.named[i,][vi]
    })
    names(mnnA) <- rownames(knnA.named)

    # get adjacency matrix
    mnn <- mnnA
    adj <- matrix(0, length(c(ctA, ctB)), length(c(ctA, ctB)))
    rownames(adj) <- colnames(adj) <- c(ctA, ctB)
    invisible(lapply(seq_len(length(mnn)), function(i) {
        # symmetric (if I'm your mutual nearest neighbor, you're mine)
        adj[names(mnn)[i],
            mnn[[i]]
            ] <<- 1
        adj[mnn[[i]],
            names(mnn)[i]
            ] <<- 1
    }))

    # reorder
    adj <- adj[rownames(pos), rownames(pos)]
    return(adj)
}


#' Nearest background neighbor
#'
#' @description Identify nearest neighbors in the background cell-type
#'
#' @param cct Vector of cell names from cell type of interest
#' @param nct Vector of cell names from background
#' @param pos Position
#' @param k Number of nearest neighbors from background for each cell from cell type of interest
#'
#' @return Boolean matrix where value = 1 if two cells are considered adjacency ie. neighbors, else 0
#'
#' @export
#'
getBnn <- function(cct, nct, pos, k) {
    knn <- RANN::nn2(pos[nct,], pos[cct,], k=k)[[1]]
    adj <- matrix(0, nrow(pos), nrow(pos))
    rownames(adj) <- colnames(adj) <- rownames(pos)
    invisible(lapply(seq_len(length(cct)), function(i) {
        adj[cct[i],nct[knn[i,]]] <<- 1
    }))
    return(adj)
}



#' Nearest neighbors across layers
#'
#' @param layers List of positions
#' @param k Number of nearest neighbors
#'
#' @return Boolean matrix where value = 1 if two cells are considered adjacency ie. neighbors, else 0
#'
#' @export
#'
getCrossLayerNeighbors <- function(layers, k=3) {
  subset <- unlist(lapply(layers, rownames))
  between <- lapply(1:(length(layers)-1), function(i) {
    pos1 <- layers[[i]]
    pos2 <- layers[[i+1]]
    ctA <- rownames(pos1)
    ctB <- rownames(pos2)
    pos <- rbind(pos1, pos2)
    if(length(ctA)==0 | length(ctB)==0) {
      wa1 <- NA
    } else {
      pos <- as.matrix(pos)
      wa1 <- getMnn(ctA, ctB, pos=pos, k=k)
    }
    return(wa1)
  })

  w <- matrix(0, nrow=length(subset), ncol=length(subset))
  rownames(w) <- colnames(w) <- subset
  ## across layers
  invisible(lapply(between, function(wa1) {
    w[rownames(wa1), colnames(wa1)] <<- wa1
  }))

  return(w)
}



#' Filter adjacency weight matrix to between two subsets of points
#'
#' @description Restrict adjacency relationships to between two subsets of points
#'
#' @param cct Cells of cell-type 1
#' @param nct Cells of cell-type 2 (assumed to be mutually exclusive with cct)
#' @param weight Adjacency weight matrix
#' @param pos Position
#' @param plot Boolean of whether to plot
#' @param ... Additional plotting parameters
#'
#' @return Boolean matrix where value = 1 if two cells are considered adjacency ie. neighbors, else 0
#'
#' @export
#'
getInterCellTypeWeight <- function(cct, nct, weight, pos=NULL, plot=FALSE, ...) {
  weightIc <- weight[c(cct, nct), c(cct, nct)]
  weightIc[nct,nct] <- 0
  weightIc[cct,cct] <- 0
  if(plot) {
    plotNetwork(pos, weightIc, ...)
    points(pos[nct,], col='blue', pch=16)
    points(pos[cct,], col='green', pch=16)
  }
  return(weightIc)
}
