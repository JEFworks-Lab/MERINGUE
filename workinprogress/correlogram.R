## sup 4

par(mfrow=c(1,1))

set.seed(0)
N <- 5^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
pos <- pos*1.5
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
pos <- jitter(pos, 1)
plotEmbedding(pos)

w <- voronoiAdjacency(pos)
rownames(w) <- colnames(w) <- rownames(pos)
plotNetwork(pos, w)
text(pos, rownames(w))

# test for a cell
i <- "cell12"

# neighbors
neighbors <- names(which(w[i,] == 1))
foo <- matrix(0, nrow(w), ncol(w))
rownames(foo) <- colnames(foo) <- rownames(w)
foo[i, neighbors] <- 1
plotNetwork(pos, foo, axes=FALSE, xlab=NA, ylab=NA, main='Neighbor Index: 1'); box()

# neighbor's neighbors
nneighbors <- names(which(colSums(w[neighbors,]==1)>=1))
nneighbors <- setdiff(nneighbors, neighbors)
foo <- matrix(0, nrow(w), ncol(w))
rownames(foo) <- colnames(foo) <- rownames(w)
foo[i, nneighbors] <- 1
plotNetwork(pos, foo, axes=FALSE, xlab=NA, ylab=NA, main='Neighbor Index: 2'); box()

# neighbor's neighbor's neighbors
nnneighbors <- names(which(colSums(w[nneighbors,]==1)>=1))
nnneighbors <- setdiff(nnneighbors, c(nneighbors, neighbors))
foo <- matrix(0, nrow(w), ncol(w))
rownames(foo) <- colnames(foo) <- rownames(w)
foo[i, nnneighbors] <- 1
plotNetwork(pos, foo, axes=FALSE, xlab=NA, ylab=NA, main='Neighbor Index: 3'); box()

####### Make into function
getNN <- function(i, K=20) {
  print(i)

  # Initial neighbors
  neighbors <- names(which(w[i,] == 1))
  foo <- matrix(0, nrow(w), ncol(w))
  rownames(foo) <- colnames(foo) <- rownames(w)
  foo[i, neighbors] <- 1
  #plotNetwork(pos, foo, main=0)
  #prevNeighbors <- c(prevNeighbors, neighbors)

  # Neighbor's neighbors
  counter <- 1
  prevNeighbors <- c(neighbors)
  nList <- list()
  nList[1] <- list(neighbors)
  while(length(neighbors)>0 & counter < K) {
    #print(counter)
    if(length(neighbors)>1) {
      nneighbors <- names(which(colSums(w[neighbors,]==1)>=1))
    } else {
      nneighbors <- names(which(w[neighbors,] == 1))
    }
    nneighbors <- setdiff(nneighbors, prevNeighbors)
    foo <- matrix(0, nrow(w), ncol(w))
    rownames(foo) <- colnames(foo) <- rownames(w)
    foo[i, nneighbors] <- 1
    #plotNetwork(pos, foo, main=counter)
    nList[counter] <- list(nneighbors)
    prevNeighbors <- c(prevNeighbors, nneighbors)
    neighbors <- nneighbors
    counter <- counter + 1
  }
  return(nList)
}

i <- "cell12"
getNN(i)



################################### Clean

getNN <- function(i, K=20) {
  print(i)

  # Initial neighbors
  neighbors <- names(which(w[i,] == 1))
  foo <- matrix(0, nrow(w), ncol(w))
  rownames(foo) <- colnames(foo) <- rownames(w)
  foo[i, neighbors] <- 1
  #plotNetwork(pos, foo, main=0)
  #prevNeighbors <- c(prevNeighbors, neighbors)

  # Neighbor's neighbors
  counter <- 1
  prevNeighbors <- c(neighbors)
  nList <- list()
  nList[1] <- list(neighbors)
  while(length(neighbors)>0 & counter < K) {
    #print(counter)
    if(length(neighbors)>1) {
      nneighbors <- names(which(colSums(w[neighbors,]==1)>=1))
    } else {
      nneighbors <- names(which(w[neighbors,] == 1))
    }
    nneighbors <- setdiff(nneighbors, prevNeighbors)
    foo <- matrix(0, nrow(w), ncol(w))
    rownames(foo) <- colnames(foo) <- rownames(w)
    foo[i, nneighbors] <- 1
    #plotNetwork(pos, foo, main=counter)
    nList[counter] <- list(nneighbors)
    prevNeighbors <- c(prevNeighbors, nneighbors)
    neighbors <- nneighbors
    counter <- counter + 1
  }
  return(nList)
}
nnWeights <- function(pos, nns, j) {
  w1 <- do.call(rbind, lapply(rownames(pos), function(i) {
    #print(i)
    cw <- rep(0, nrow(pos))
    names(cw) <- rownames(pos)
    if(length(nns[[i]])>=j) {
      cnbs <- nns[[i]][[j]]
      if(!identical(cnbs, character(0))) {
        cw[cnbs] <- 1
      }
    }
    return(cw)
  }))
  return(w1)
}


par(mfrow=c(1,1))

set.seed(0)
N <- 20^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
pos <- pos*1.5
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)

pos[,1] <- expPos(pos[,1])
pos[,2] <- expPos(pos[,2])
plot(pos)

pos <- jitter(pos, 1)
plotEmbedding(pos)

w <- voronoiAdjacency(pos)
rownames(w) <- colnames(w) <- rownames(pos)
plotNetwork(pos, w)
#text(pos, rownames(w))

# neighbors neighbors for each cell
nns <- lapply(rownames(pos), getNN)
names(nns) <- rownames(pos)


############## Sample gene expression
par(mfrow=c(3,2))

gexp <- rep(-1, nrow(pos))
vi <- pos[,1]>0.6*max(pos[,1]) & pos[,2]>0.6*max(pos[,2])
gexp[vi] <- 1
vi2 <- pos[,1]<0.3*max(pos[,1]) & pos[,2]<0.3*max(pos[,2])
gexp[vi2] <- 1
names(gexp) <- rownames(pos)
plotEmbedding(pos, col=gexp, cex=1, main='Simulated Gene A')
x <- moranTest(gexp, w)$observed
x <- c(x, unlist(lapply(1:10, function(j) {
  w1 <- nnWeights(pos, nns, j)
  moranTest(gexp, w1)$observed
})))

slisa <- lisaTest(gexp, w)$observed*sign(gexp)
names(slisa) <- names(gexp)
slisa <- scale(slisa)[,1]
slisa[slisa < -0.5] <- -0.5
slisa[slisa > 0.5] <- 0.5
plotEmbedding(pos, col=slisa, main='sLISA', cex=1)


#### spatial correlogram

gexp <- scale(pos[,1]*pos[,2])[,1]
names(gexp) <- rownames(pos)
plotEmbedding(pos, col=gexp, cex=1, main='Simulated Gene B')
x <- moranTest(gexp, w)$observed
x <- c(x, unlist(lapply(1:10, function(j) {
  w1 <- nnWeights(pos, nns, j)
  moranTest(gexp, w1)$observed
})))
x1 <- x
x1

plot(1:length(x1), x1, type="l", col='blue', xlab='Neighbor Index', ylab='Observed Moran\'s Statistic', main='Correlogram')


## alternating pattern?
gexp <- rep(-1, nrow(pos))
vi <- seq(1, 20^2, 40)
vii <- lapply(2:length(vi), function(i) seq(vi[i-1], vi[i]-1))
viii <- unlist(vii[seq(1, length(vi), 2)])
gexp[viii] <- 1
names(gexp) <- rownames(pos)
plotEmbedding(pos, col=gexp, cex=1, main='Simulated Gene C')

moranTest(gexp, w, alternative='two.sided')
x <- moranTest(gexp, w)$observed
x <- c(x, unlist(lapply(1:10, function(j) {
  w1 <- nnWeights(pos, nns, j)
  moranTest(gexp, w1)$observed
})))
x3 <- x
plot(1:length(x3), x3, type="l", col='blue', xlab='Neighbor Index', ylab='Observed Moran\'s Statistic', main='Correlogram')


