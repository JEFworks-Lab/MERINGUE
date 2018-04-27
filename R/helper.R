

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


#' Moran's I compute from scratch
#' x is value
#' w is adjacency matrix (weights)
moranTest <- function(x, w) {
  ## Compute from scratch
  #w <- w / sum(w)
  #n <- length(x)
  #z <- as.vector((x - mean(x)) / sd(x))
  #as.vector(z %*% w %*% (z * sqrt(n / (n-1))))

  ## Just use Ape's package since includes empirical p-values
  ## Only look for greater autocorrelation / increased clustering
  unlist(ape::Moran.I(as.vector(x), w, alternative='greater'))
}

#' Permutation test to assess for significance of Moran's I
moranPermutationTest <- function(z, w, N=1e4, seed=0, ncores=parallel::detectCores()-1, plot=FALSE, ...) {
  # Set seed for reproducibility
  set.seed(seed)
  # Compute Moran's I
  stat <- moranTest(z, w)['observed']
  # Simulate null distribution
  sim <- unlist(parallel::mclapply(seq_len(N), function(i) {
    moranTest(sample(z, length(z), replace=TRUE), w)['observed']
  }, mc.cores=ncores))
  p.value <- mean((all <- c(stat, sim)) >= stat)
  if(plot) {
    hist(sim, sub=paste("p =", round(p.value, 4)), xlim=range(all), ...)
    abline(v = stat, col="red", lty=3, lwd=2)
  }
  results <- unlist(data.frame('observed'=stat, 'N'=N, 'p.value'=p.value))
  return(results)
}


#' Gridded bivariate interpolation
interpolate <- function(pos, gexp, binSize=100, col=colorRampPalette(c('blue', 'white', 'red'))(100), plot=TRUE) {
  z <- gexp
  x <- pos[,1]
  y <- pos[,2]
  int <- akima::interp(x, y, z, nx=binSize, ny=binSize, linear=FALSE)
  if(plot) {
    plot(pos, col=map2col(z), pch=16, cex=2, axes=FALSE, frame.plot=TRUE, xlab=NA, ylab=NA)
    image(int, col=col, axes=FALSE, frame.plot=TRUE)
  }
  return(int)
}


#' 2D kernel density estimation
densityPlot <- function(pos) {
  x <- pos[,1]
  y <- pos[,2]
  dens <- MASS::kde2d(x, y)
  persp(dens, phi = 30, theta = 20, d = 5)
}
