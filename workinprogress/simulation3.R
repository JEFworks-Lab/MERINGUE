## Spatial cross correlation simulation
library(MERingue)

set.seed(0)
N <- 100
pos <- cbind(rnorm(N), rnorm(N))
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')

weight <- voronoiAdjacency(pos, filterDist=2)
## Simulate spatially cross correlated but not necessarily correlated expression due to different cell-types being affected
ctA <- sample(rownames(pos), N/2)
ctB <- setdiff(rownames(pos), ctA)

par(mfrow=c(1,1))
plotNetwork(pos, weight, line.col='grey')
points(pos[ctA,], col='red')
points(pos[ctB,], col='blue')

gexpA <- pos[,2]
gexpA[ctB] <- 0
gexpB <- -pos[,2]
gexpB[ctA] <- 0
#gexpA[gexpA < 0] <- 0
#gexpB[gexpB < 0] <- 0

plotEmbedding(pos, colors=gexpA, main='gene A', cex=1)
plotEmbedding(pos, colors=gexpB, main='gene B', cex=1)

plot(gexpA, gexpB) # no correlation
ct <- cor.test(gexpA, gexpB)
text(paste0('R^2 ', ct$estimate), col='red', x=-1, y=1)
text(paste0('p-value ', ct$p.value), col='red', x=-1, y=1.2)

## spatial cross correlation
I <- spatialIntraCrossCor(gexpA, gexpB, weight)
spatialCrossCorMatrix(rbind(gexpA, gexpB), weight)
I

# random permutation
perm <- sapply(1:1000, function(i) spatialIntraCrossCor(sample(gexpA), sample(gexpB), weight))
hist(perm, breaks=50, xlim=c(-1,1))
abline(v=I, col='red')
pv <- 2 * pnorm(abs(I), mean(perm), sd(perm), lower.tail=FALSE) # two-sided p-value
text(paste0('p-value ', pv), col='red', x=0, y=50)

plotNeighborBoxplot(gexpA, gexpB, ctA, ctB, weight)

# toroidal shift model
perm <- sapply(1:1000, function(i) {
  ppos <- rtorShift(pos, seed = i, k=1)
  pweight <- voronoiAdjacency(ppos, filterDist=1)

  #par(mfrow=c(2,2))
  #plotNetwork(ppos, pweight, line.col='grey')
  #points(ppos[ctA,], col='red')
  #points(ppos[ctB,], col='blue')
  #plotEmbedding(ppos, colors=gexpA, main='gene A', cex=1)
  #plotEmbedding(ppos, colors=gexpB, main='gene B', cex=1)

  spatialIntraCrossCor(gexpA, gexpB, pweight)
})
par(mfrow=c(1,1))
hist(perm, breaks=50, xlim=c(mean(perm)-0.5,mean(perm)+0.5))
abline(v=I, col='red')
I
pv <- 2 * pnorm(abs(I), mean(perm), sd(perm), lower.tail=FALSE) # two-sided p-value
text(paste0('p-value ', pv), col='red', x=mean(perm), y=50)

##





############# Test on OB
set.seed(0)
library(MERingue)
data(mOB)

# Clean and normalize data
pos <- mOB$pos
cd <- mOB$counts
counts <- cleanCounts(cd, min.reads=10, min.lib.size=10)
pos <- pos[colnames(counts),]
mat <- normalizeCounts(counts, log=FALSE)

w <- voronoiAdjacency(pos, njitter = 10, ajitter = 5, filterDist = 1.75)
par(mfrow=c(1,1))
plotNetwork(pos, w)
I <- getSpatialPatterns(mat, w)
results.filter <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
length(results.filter)

results.sig <- results.filter
results.sig
lapply(results.sig, function(g1) {
  results.g1 <- lapply(results.sig, function(g2) {
    print(g2)

    gexpA <- mat[g1,]
    gexpB <- mat[g2,]

    #plot(gexpA, gexpB)
    #plotEmbedding(pos, colors=scale(gexpA)[,1], main=g1, cex=2)
    #plotEmbedding(pos, colors=scale(gexpB)[,1], main=g2, cex=2)

    # regular correlation test
    ct <- cor.test(gexpA, gexpB)
    ct

    # spatiall correlation test
    SCI <- spatialIntraCrossCor(gexpA, gexpB, w)
    perm <- sapply(1:100, function(i) spatialIntraCrossCor(sample(gexpA), sample(gexpB), w))

    #hist(perm, breaks=50, xlim=c(-1,1))
    #abline(v=SCI, col='red')

    SCIpv <- 2 * pnorm(abs(SCI), mean(perm), sd(perm), lower.tail=FALSE) # two-sided p-value

    data.frame('cor'=ct$estimate, 'cor.pv'=ct$p.value, 'SCI'=SCI, 'SCI.pv'=SCIpv)
  })
  names(results.g1) <- results.sig
  do.call(rbind, results.g1)
})

