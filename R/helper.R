

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
plotNetwork <- function(pos, adj, col='black', line.col='red', line.power=1, ...) {
  plot(pos, pch=16, col=col)
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
#' weight is adjacency matrix (weights)
#' use log.p since many p-values very small (avoid 0s)
moranTest <- function (x, weight, na.rm = FALSE, alternative = "greater") {

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- length(x)

  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

  nas <- is.na(x)
  if (any(nas)) {
    if (na.rm) {
      x <- x[!nas]
      N <- length(x)
      weight <- weight[!nas, !nas]
    }
    else {
      stop("'x' has missing values")
    }
  }

  # first moment
  ei <- -1/(N - 1)

  # unitization
  rs <- rowSums(weight)
  rs[rs == 0] <- 1
  weight <- weight/rs

  # Moran's I
  W <- sum(weight)
  z <- x - mean(x)
  cv <- sum(weight * z %o% z)
  v <- sum(z^2)
  obs <- (N/W) * (cv/v)

  # second moment
  W.sq <- W^2
  N.sq <- N^2
  S1 <- 0.5 * sum((weight + t(weight))^2)
  S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
  S3 <- (sum(z^4)/N)/(v/N)^2
  S4 <- (N.sq - 3*N + 3)*S1 - N*S2 + 3*W.sq
  S5 <- (N.sq - N)*S1 - 2*N*S2 + 6*W.sq
  ei2 <- (N*S4 - S3*S5)/((N - 1)*(N - 2)*(N - 3) * W.sq)

  # standard deviation
  sdi <- sqrt(ei2 - (ei)^2)

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  pv <- pnorm(obs, mean = ei, sd = sdi)
  if (alternative == "two.sided") {
    if (obs <= ei) {
      pv <- 2 * pv
    } else {
      pv <- 2 * (1 - pv)
    }
  }
  if (alternative == "greater") {
    pv <- 1 - pv
  }

  return(list(observed = obs, expected = ei, sd = sdi, p.value = pv))
}


moranSimple <- function(x, weight, na.rm = FALSE) {
  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- length(x)

  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

  nas <- is.na(x)
  if (any(nas)) {
    if (na.rm) {
      x <- x[!nas]
      N <- length(x)
      weight <- weight[!nas, !nas]
    }
    else {
      stop("'x' has missing values")
    }
  }

  # scale weights
  rs <- rowSums(weight)
  rs[rs == 0] <- 1
  weight <- weight/rs

  # Moran's I
  W <- sum(weight)
  m <- mean(x)
  y <- x - m
  cv <- sum(weight * y %o% y)
  v <- sum(y^2)
  obs <- (N/W) * (cv/v)

  return(obs)
}


#' Permutation test to assess for significance of Moran's I
moranPermutationTest <- function(z, w, na.rm = FALSE, alternative = "two.sided", N=1e4, seed=0, ncores=parallel::detectCores()-1, plot=FALSE, ...) {
  # Set seed for reproducibility
  set.seed(seed)
  # Compute Moran's I
  stat <- moranSimple(z, w)
  # Simulate null distribution
  sim <- unlist(parallel::mclapply(seq_len(N), function(i) {
    moranSimple(sample(z, length(z), replace=TRUE), w)
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
interpolate <- function(pos, gexp, binSize=100, cex=1, col=colorRampPalette(c('blue', 'white', 'red'))(100), plot=TRUE, ...) {
  z <- gexp
  x <- pos[,1]
  y <- pos[,2]
  int <- akima::interp(x, y, z, nx=binSize, ny=binSize, linear=FALSE)
  if(plot) {
    plot(pos, col=map2col(z), pch=16, cex=cex, axes=FALSE, frame.plot=TRUE, xlab=NA, ylab=NA, ...)
    image(int, col=col, axes=FALSE, frame.plot=TRUE, ...)
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




#' Winsorize expression values to prevent outliers
winsorize <- function (x, fraction=.05) {
   if(length(fraction) != 1 || fraction < 0 ||
         fraction > 0.5) {
      stop("bad value for 'fraction'")
   }
   lim <- quantile(x, probs=c(fraction, 1-fraction))
   x[ x < lim[1] ] <- lim[1]
   x[ x > lim[2] ] <- lim[2]
   x
}
