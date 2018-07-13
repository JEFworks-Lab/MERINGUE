#' Moran's I compute from scratch
#' x is value
#' weight is adjacency matrix (weights)
#' use log.p since many p-values very small (avoid 0s)
#' DEPRECATED -> use moranTest_C
#'
#' @export
#'
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


#' Permutation test to assess for significance of Moran's I
#' @export
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



#' Winsorize expression values to prevent outliers
#' @export
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




#' inter-cell-type cross cor
#' @export
spatialCrossCor <- function(gexpA, gexpB, groupA, groupB, weight=NULL, pos=NULL, k=3) {
    # make ctA the smaller group
    if(length(groupA) < length(groupB)) {
        ctA <- groupA
        ctB <- groupB
    } else {
        ctA <- groupB
        ctB <- groupA
    }

    # restrict to expression of gene A in group A
    # and expression of gene B in group B
    x <- gexpA
    y <- gexpB
    x[groupB] <- NA
    y[groupA] <- NA

    if(is.null(weight)) {
        weight <- getMnn(ctA, ctB, pos, k)
    }

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
    SCI <- (N/W) * (cv/v)

    return(SCI)
}


#' intra-cell-type cross cor
#' @export
spatialIntraCrossCor <- function(x, y, weight) {
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
    SCI <- (N/W) * (cv/v)
    SCI
}


# Local indicators of spatial association
#' @export
lisaTest <- function (x, weight, na.rm=FALSE, alternative = "greater") {
    # unitization
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight/rs

    # calculate Ii
    n <- length(x)
    xx <- mean(x, na.rm = na.rm)
    z <- x - xx
    s2 <- sum(z^2, na.rm = na.rm)/n
    lz <- apply(weight, 1, function(w) sum(w*z))
    Ii <- (z/s2) * lz

    Wi <- rowSums(weight)
    E.Ii <- -Wi/(n - 1)

    b2 <- (sum(z^4, na.rm = na.rm)/n)/(s2^2)
    Wi2 <- apply(weight, 1, function(w) sum(w^2))
    A <- (n - b2)/(n - 1)
    B <- (2 * b2 - n)/((n - 1) * (n - 2))
    Sd.Ii <- sqrt(A * Wi2 + B * (Wi^2 - Wi2) - E.Ii^2)

    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    pv <- pnorm(Ii, mean = E.Ii, sd = Sd.Ii)
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

    return(data.frame(observed = Ii, expected = E.Ii, sd = Sd.Ii, p.value = pv))
}


#' Calculates a spatial cross correlation matrix
#'
#' @export
spatialCrossCorMatrix <- function(sigMat, w) {
  scor <- spatialCrossCorMatrix_C(as.matrix(sigMat), w)
  colnames(scor) <- rownames(scor) <- rownames(sigMat)
  return(scor)
}
