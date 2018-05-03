# calculate moran's I statistic from scratch

x = value
w = adj
ape::Moran.I(x, w, alternative='greater')

moran <- function (x, weight, na.rm = FALSE, alternative = "two.sided") {

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

  # second moment
  W.sq <- W^2
  N.sq <- N^2
  S1 <- 0.5 * sum((weight + t(weight))^2)
  S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
  S3 <- (sum(y^4)/N)/(v/N)^2
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
