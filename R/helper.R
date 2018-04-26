#' Helper function to map values to colors
#' Source: https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
map2col <- function(x, pal=colorRampPalette(c('blue', 'white', 'red'))(100), limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

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

#' Plot neighbor network
#' https://stackoverflow.com/questions/43879347/plotting-a-adjacency-matrix-using-pure-r
plotNetwork <- function(pos, adj, line.col='red', line.power=1, ...) {
  plot(pos, pch=16)
  idx <- which(adj>0, arr.ind = T)
  for(i in seq_len(nrow(idx))) {
    lines(
      c(pos[idx[i,1],1], pos[idx[i,2],1]),
      c(pos[idx[i,1],2], pos[idx[i,2],2]),
      col=line.col,
      lwd=adj[idx]^line.power
      )
  }
}
