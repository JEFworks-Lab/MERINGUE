#' Helper function to get p-value from value and distribution mean and standard deviation
#'
#' @param obs Observed value
#' @param ei Mean (Expected value)
#' @param sdi Standard deviation
#' @param alternative "two.sided", "less", or "greater"
#'
#' @return P-value
#'
#' @examples
#' # Probability of observing 5 while expecting 0 with a standard deviation of 1
#' getPv(5, 0, 1, 'greater')
#'
#' @export
#'
getPv <- function(obs, ei, sdi, alternative) {
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
  if (alternative == "lesser") {
    pv <- pv
  }
  return(pv)
}


#' Spatialy autocorrelation by Moran's I in C++
#'
#' @param x Feature value
#' @param weight Adjacency weight matrix
#' @param alternative "two.sided", "less", or "greater"
#'
#' @return Observed Moran's I statistic, the expected statistic under the null hypothesis of no spatial autocorrelation, the standard deviation under the null hypothesis, and p-value
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' weight <- voronoiAdjacency(pos)
#' gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)['Camk4',]
#' plotEmbedding(pos, colors=scale(gexp)[,1], cex=3)
#' moranTest(gexp, weight)
#'
#' @export
#'
moranTest <- function (x, weight, alternative = "greater") {

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- length(x)
  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

  if(sum(rownames(weight) %in% names(x)) != nrow(weight)) {
    stop('Names in feature vector and adjacency weight matrix do not agree.')
  }
  x <- x[rownames(weight)]

  output <- moranTest_C(x, weight)
  output <- output[,1]
  names(output) <- c('observed', 'expected', 'sd')
  output['p.value'] <- getPv(output[1], output[2], output[3], alternative)
  return(output)
}


#' Spatialy autocorrelation by Moran's I
#' DEPRECATED! Use moranTest()
#'
#' @param x Feature value
#' @param weight Adjacency weight matrix
#' @param alternative "two.sided", "less", or "greater"
#'
#' @return Observed Moran's I statistic, the expected statistic under the null hypothesis of no spatial autocorrelation, the standard deviation under the null hypothesis, and p-value
#'
#' @export
#'
moranTest_DEPRECATED <- function (x, weight, alternative = "greater") {

  warning('DEPRECATED! Use moranTest()')

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- length(x)
  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

  if(sum(rownames(weight) %in% names(x)) != nrow(weight)) {
    stop('Names in feature vector and adjacency weight matrix do not agree.')
  }
  x <- x[rownames(weight)]

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

  # p-value
  pv <- getPv(obs, ei, sdi, alternative)

  return(c(observed = obs, expected = ei, sd = sdi, p.value=pv))
}


#' Derive Moran's I p-value using permutation testing
#'
#' @param z Feature value
#' @param w Adjacency weight matrix
#' @param alternative "two.sided", "less", or "greater"
#' @param N Number of permutations
#' @param seed Random seed
#' @param ncores Number of cores for parallel processing
#' @param plot Plot permutated distribution
#' @param ... Additional parameters to pass to histogram plotting
#'
#' @return Observed Moran's I statistic, the expected statistic under the null hypothesis of no spatial autocorrelation, the standard deviation under the null hypothesis, permutation p-value, and number of permutations
#'
#' @export
#'
moranPermutationTest <- function(z, w, alternative = "greater", N=1e4, seed=0, ncores=1, plot=FALSE, ...) {

  # Set seed for reproducibility
  set.seed(seed)

  # Compute Moran's I
  stat <- moranSimple(z, w)

  # Simulate null distribution
  sim <- unlist(parallel::mclapply(seq_len(N), function(i) {
    moranSimple(sample(z, length(z), replace=TRUE), w)
  }, mc.cores=ncores))
  sim[is.nan(sim)] <- NA
  sim <- na.omit(sim)
  all <- c(stat, sim)

  # Compute permutation p-value
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if(alternative == 'greater') {
    p.value <- mean(all >= stat)
  }
  if(alternative == 'less') {
    p.value <- mean(all <= stat)
  }
  if(alternative == 'two.sided') {
    p.value <- mean(abs(all) >= abs(stat))
  }

  # Plot null distribution
  if(plot) {
    hist(sim, sub=paste("p =", round(p.value, 4)), xlim=range(all), ...)
    abline(v = stat, col="red", lty=3, lwd=2)
  }

  results <- c('observed'=stat, 'expected'=mean(sim), 'sd'=sd(sim), 'p.value'=p.value, 'N'=N)
  return(results)
}
moranSimple <- function(x, weight) {
  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- length(x)
  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

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

  return(obs)
}


#' Calculates a spatial cross correlation for all pairs
#'
#' @param mat Matrix of feature values
#' @param weight Adjacency weight matrix
#'
#' @return Matrix of spatial cross correlationf or all pairs
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' weight <- voronoiAdjacency(pos)
#' gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)
#' sccMat <- spatialCrossCorMatrix(gexp[1:5,], weight)
#'
#' @export
#'
spatialCrossCorMatrix <- function(mat, weight) {

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- ncol(mat)
  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

  if(sum(rownames(weight) %in% colnames(mat)) != nrow(weight)) {
    stop('Names in feature vector and adjacency weight matrix do not agree.')
  }
  mat <- mat[, rownames(weight)]

  scor <- spatialCrossCorMatrix_C(as.matrix(mat), weight)
  colnames(scor) <- rownames(scor) <- rownames(mat)
  return(scor)
}


#' Calculates spatial cross correlation for two features
#'
#' @param x Feature 1 value
#' @param y Feature 2 value
#' @param weight Adjacency weight matrix
#'
#' @return Spatial cross correlation statistic
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' weight <- voronoiAdjacency(pos)
#' gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)
#' scc <- spatialCrossCor(gexp[1,], gexp[2,], weight)
#'
#' @export
#'
spatialCrossCor <- function(x, y, weight) {
  if(length(x) != length(y)) {
    stop("'x', and 'y' must be equal in length")
  }
  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }
  N <- length(x)
  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x' and 'y'")
  }

  if(sum(rownames(weight) %in% names(x)) != nrow(weight)) {
    stop('Names in feature vector and adjacency weight matrix do not agree.')
  }
  x <- x[rownames(weight)]
  y <- y[rownames(weight)]

  # scale weights
  rs <- rowSums(weight)
  rs[rs == 0] <- 1
  weight <- weight/rs

  # compute spatial cross correlation
  N <- length(x)
  W <- sum(weight)
  dx <- x - mean(x, na.rm=TRUE)
  dy <- y - mean(y, na.rm=TRUE)

  cv1 <- dx %o% dy
  cv2 <- dy %o% dx
  cv1[is.na(cv1)] <- 0
  cv2[is.na(cv2)] <- 0

  cv <- sum(weight * ( cv1 + cv2 ), na.rm=TRUE)
  v <- sqrt(sum(dx^2, na.rm=TRUE) * sum(dy^2, na.rm=TRUE))
  SCC <- (N/W) * (cv/v) / 2.0

  return(SCC)
}


#' Tests for inter-cell-type spatial cross-correlation between gene A in group A and gene B in group B
#'
#' @param gexpA Expression for gene A
#' @param gexpB Expression for gene B
#' @param groupA Cells in group A
#' @param groupB Cells in group B
#' @param weight Adjacency weight matrix
#'
#' @return statistic of how well gene expression A in cells of group A are correlated in space with gene expression B in cells of group B
#'
#' @examples
#' # Simulate data
#' set.seed(0)
#' N <- 100
#' pos <- cbind(rnorm(N), rnorm(N))
#' rownames(pos) <- paste0('cell', 1:N)
#' colnames(pos) <- c('x', 'y')
#' weight <- voronoiAdjacency(pos)
#' ctA <- sample(rownames(pos), N/2)
#' ctB <- setdiff(rownames(pos), ctA)
#' par(mfrow=c(1,1))
#' plotNetwork(pos, weight, line.col='grey')
#' points(pos[ctA,], col='red')
#' points(pos[ctB,], col='blue')
#' gexpA <- pos[,2]
#' gexpA[ctB] <- 0
#' gexpB <- -pos[,2]
#' gexpB[ctA] <- 0
#' plotEmbedding(pos, col=gexpA)
#' plotEmbedding(pos, col=gexpB)
#' interCellTypeSpatialCrossCor(gexpA, gexpB, ctA, ctB, weight)
#'
#' @export
#'
interCellTypeSpatialCrossCor <- function(gexpA, gexpB, groupA, groupB, weight) {
    # restrict to expression of gene A in group A
    # and expression of gene B in group B
    x <- gexpA
    y <- gexpB
    x[groupB] <- NA
    y[groupA] <- NA

    # scale weights
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight/rs

    # compute spatial cross correlation
    N <- length(x)
    W <- sum(weight)
    dx <- x - mean(x, na.rm=TRUE)
    dy <- y - mean(y, na.rm=TRUE)

    cv1 <- dx %o% dy
    cv2 <- dy %o% dx
    cv1[is.na(cv1)] <- 0
    cv2[is.na(cv2)] <- 0

    cv <- sum(weight * ( cv1 + cv2 ), na.rm=TRUE)
    v <- sqrt(sum(dx^2, na.rm=TRUE) * sum(dy^2, na.rm=TRUE))
    SCI <- (N/W) * (cv/v) / 2.0

    return(SCI)
}


#' Local indicators of spatial association
#'
#' @param x Feature value
#' @param weight Adjacency weight matrix
#' @param alternative "two.sided", "less", or "greater"
#'
#' @return Data frame of observed, expected, standard deviation, and p-value for each point
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' weight <- voronoiAdjacency(pos)
#' gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)['Camk4',]
#' lisa <- lisaTest(gexp, weight)
#'
#' @export
#'
lisaTest <- function(x, weight, alternative = "greater") {

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix")
  }

  N <- length(x)
  if (nrow(weight) != N) {
    stop("'weight' must have as many rows as observations in 'x'")
  }

  if(sum(rownames(weight) %in% names(x)) != nrow(weight)) {
    stop('Names in feature vector and adjacency weight matrix do not agree.')
  }
  x <- x[rownames(weight)]

    # unitization
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight/rs

    # calculate Ii
    n <- length(x)
    xx <- mean(x, na.rm = TRUE)
    z <- x - xx
    s2 <- sum(z^2, na.rm = TRUE)/n
    lz <- apply(weight, 1, function(w) sum(w*z))
    Ii <- (z/s2) * lz

    Wi <- rowSums(weight)
    E.Ii <- -Wi/(n - 1)

    b2 <- (sum(z^4, na.rm = TRUE)/n)/(s2^2)
    Wi2 <- apply(weight, 1, function(w) sum(w^2))
    A <- (n - b2)/(n - 1)
    B <- (2 * b2 - n)/((n - 1) * (n - 2))
    Sd.Ii <- sqrt(A * Wi2 + B * (Wi^2 - Wi2) - E.Ii^2)

    # p-value
    pv <- getPv(Ii, E.Ii, Sd.Ii, alternative)

    return(data.frame(observed = Ii, expected = E.Ii, sd = Sd.Ii, p.value = pv))
}


#' Signed local indicators of spatial association to identify regions driving global autocorrelation
#'
#' @param x Feature value
#' @param weight Adjacency weight matrix
#' @param alternative "two.sided", "less", or "greater"
#'
#' @return Signed LISA statistic for each point
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' weight <- voronoiAdjacency(pos)
#' gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)['Camk4',]
#' plotEmbedding(pos, colors=scale(gexp)[,1], cex=3)
#' slisa <- signedLisa(gexp, weight)
#' plotEmbedding(pos, colors=slisa,
#'    gradientPalette=colorRampPalette(c('darkgreen', 'white', 'darkorange'))(100))
#'
#' @export
#'
signedLisa <- function(x, weight, alternative = "greater") {
  lisa <- lisaTest(x, weight, alternative)
  slisa <- sign(scale(x[rownames(lisa)], center=TRUE)[,1])*lisa$observed
  slisa <- slisa[rownames(weight)]
  return(slisa)
}

