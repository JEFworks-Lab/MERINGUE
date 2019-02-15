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


#' Adjacency weight matrix by voronoi adjacency
#'
#' @description Adjacency weight matrix by voronoi adjacency. Modified from Barry Rowlingson and Alison Hale's voronoi_adjacency function in the caramellar package
#'
#' @param pos Position
#' @param filterDist Euclidean distance beyond which two cells cannot be considered neighbors
#' @param nDummy Number of dummy points to use for edge cases
#' @param plot Boolean of whether to plot
#'
#' @return Boolean matrix where value = 1 if two cells are considered adjacency ie. neighbors, else 0
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' w <- voronoiAdjacency(pos, plot=TRUE)
#'
#' @export
#'
voronoiAdjacency <- function(pos, filterDist = NA, nDummy = 3, plot=FALSE){

  if(ncol(pos)>2) {
    stop('2D tesselation only')
  }
  if(is.null(colnames(pos))) {
    colnames(pos) <- c('x', 'y')
  }
  if(is.null(rownames(pos))) {
    rownames(pos) <- seq_len(nrow(pos))
  }

  ## deldir doesn't handle edge cases well
  pos0 <- pos
  ## create dummy edges
  helper <- function(pos, nDummy) {
    t <- max((max(pos[,1])-min(pos[,1]))/nDummy, (max(pos[,2])-min(pos[,2]))/nDummy)
    x1 <- min(pos[,1])-t
    x2 <- max(pos[,1])+t
    y1 <- min(pos[,2])-t
    y2 <- max(pos[,2])+t
    dummy <- rbind(
      cbind(seq(x1, x2, by=(x2-x1)/nDummy), y1),
      cbind(seq(x1, x2, by=(x2-x1)/nDummy), y2),
      cbind(x1, seq(y1, y2, by=(y2-y1)/nDummy)),
      cbind(x2, seq(y1, y2, by=(y2-y1)/nDummy))
    )
    colnames(dummy) <- c('x', 'y')
    rownames(dummy) <- paste0('dummy', 1:nrow(dummy))
    pos <- rbind(pos, dummy)

    ## find adjacencies
    makeixy <- function(data, formula, scale){
      m = model.frame(formula, data=data)
      if(ncol(m)!=3){
        stop("incorrect adjacency formula: id~x+y needed")
      }
      names(m)=c("id","x","y")
      m[,2]=m[,2]/scale
      m[,3]=m[,3]/scale
      m
    }
    data <- as.data.frame(cbind('i'=seq_len(nrow(pos)), pos))
    data <- makeixy(data, formula=i~x+y, scale=1)
    return(list('pos'=pos, 'data'=data))
  }
  dd1 <- NULL
  i <- nDummy
  while(is.null(dd1)) {
    warning(paste0('Deldir with nDummy ', i+1, '.'))
    data.info <- helper(pos0, i)
    data <- data.info$data
    pos <- data.info$pos
    dd1 <- deldir::deldir(data$x,
                          data$y,
                          suppressMsge=TRUE, sort=FALSE,
                          plotit = plot, col='grey')
    i = i+1
  }
  if(plot) { box() }

  ## final
  data.info <- helper(pos, i)
  data <- data.info$data
  pos <- data.info$pos

  ## create distance matrix
  P <- nrow(data) # number of rows
  D1 <- matrix(0,P,P)
  rownames(D1) <- rownames(pos)
  colnames(D1) <- rownames(pos)
  D1[as.matrix(dd1$delsgs[,c("ind1","ind2")])] = sqrt((dd1$delsgs[,c("x1")]-dd1$delsgs[,c("x2")])^2+(dd1$delsgs[,c("y1")]-dd1$delsgs[,c("y2")])^2)
  D1[as.matrix(dd1$delsgs[,c("ind2","ind1")])] = sqrt((dd1$delsgs[,c("x1")]-dd1$delsgs[,c("x2")])^2+(dd1$delsgs[,c("y1")]-dd1$delsgs[,c("y2")])^2)
  D <- D1

  ## filter by distance
  if(!is.na(filterDist)) {
    D[D>filterDist] = 0
  }
  ## filter to original set of points
  D <- D[rownames(pos0), rownames(pos0)]

  ## binarize
  D[D>0] <- 1

  ## plot
  if(plot) {
    idx <- which(D>0, arr.ind = T)
    for(i in seq_len(nrow(idx))) {
      lines(
        c(pos0[idx[i,1],1], pos0[idx[i,2],1]),
        c(pos0[idx[i,1],2], pos0[idx[i,2],2]),
        col='red'
      )
    }
    points(pos0, pch=16)
  }


  return(D);
}
