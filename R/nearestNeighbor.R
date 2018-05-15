#' K nearest neighbors
getAdj <- function(mat, k) {
  ## nearest neighbors include self so add 1
  knn <- RANN::nn2(mat, k=k+1)[[1]]
  knn <- knn[, -1]
  ## convert to adjacency matrix
  adj <- matrix(0, nrow(mat), nrow(mat))
  rownames(adj) <- colnames(adj) <- rownames(mat)
  invisible(lapply(seq_len(nrow(mat)), function(i) {
    adj[i,rownames(mat)[knn[i,]]] <<- 1
  }))
  return(adj)
}

#' Mutual nearest neighbors
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
getBnn <- function(cct, nct, pos, k) {
    knn <- RANN::nn2(pos[nct,], pos[cct,], k=k)[[1]]
    adj <- matrix(0, nrow(pos), nrow(pos))
    rownames(adj) <- colnames(adj) <- rownames(pos)
    invisible(lapply(seq_len(length(cct)), function(i) {
        adj[cct[i],nct[knn[i,]]] <<- 1
    }))
    return(adj)
}




