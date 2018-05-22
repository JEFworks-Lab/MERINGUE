## Spatial cross correlation simulation

N <- 100
pos <- cbind(rnorm(N), rnorm(N))
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')

weight <- getAdj(pos, k=3)
plotNetwork(pos, weight, line.col='grey')

## Simulate spatially cross correlated but not necessarily correlated expression due to different cell-types being affected
ctA <- sample(rownames(pos), N/2)
ctB <- setdiff(rownames(pos), ctA)
points(pos[ctA,], col='red')
points(pos[ctB,], col='blue')

gexpA <- pos[,2]
gexpA[ctB] <- 0
gexpB <- -pos[,2]
gexpB[ctA] <- 0

plotEmbedding(pos, colors=gexpA, main='gene A')
plotEmbedding(pos, colors=gexpB, main='gene B')

plot(gexpA, gexpB) # no correlation

## spatial cross correlation
I <- spatialIntraCrossCor(gexpA, gexpB, weight)
I

# random permutation
perm <- sapply(1:1000, function(i) spatialIntraCrossCor(sample(gexpA), sample(gexpB), weight))
hist(perm, breaks=50, xlim=c(-1,1))
abline(v=I, col='red')

2 * pnorm(abs(I), mean(perm), sd(perm), lower.tail=FALSE) # two-sided p-value


############# Test on OB
library(MERingue)

# Data
data(cd)
data(pos)

cd <- as.matrix(cd[rownames(pos),])
sample_info <- pos
pos <- pos[,1:2]
hist(cd)
counts <- cleanCounts(t(as.matrix(cd)), min.reads=100, min.detected=100)
mat <- normalizeCounts(counts)

w <- getSpatialWeights(pos, klist=3)
par(mfrow=c(1,1))
plotNetwork(pos, w)
I <- getSpatialPatterns(mat, w)

results <- I
vi <- results$p.value < 0.05
vi[is.na(vi)] <- FALSE
table(vi)
results.sig <- rownames(results)[vi]

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

