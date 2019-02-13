## Double check that Moran's I p-value calculated from moments is the same as permutation

library(MERingue)
library(Matrix)

## load mOB data
data(mOB)
# Clean and normalize data
pos <- mOB$pos
cd <- mOB$counts

counts <- cleanCounts(cd, min.reads=10, min.lib.size=10)
pos <- pos[colnames(counts),]
mat <- normalizeCounts(counts, log=FALSE)

w <- voronoiAdjacency(pos, njitter = 10, ajitter = 5, filterDist = 1.75)

## check for one gene
i <- 1
moranTest(mat[i,], w)$p.value
moranPermutationTest(mat[i,], w, plot=TRUE)$p.value

## do for N genes
set.seed(0)
N <- 100
genes <- sample(rownames(mat), N)

p.comp <- do.call(rbind, lapply(genes, function(i) {
  print(i)
  p1 <- moranTest(mat[i,], w)$p.value
  p2 <- moranPermutationTest(mat[i,], w, plot=FALSE)$p.value
  return(c(p1, p2))
}))

## assess correspondence
par(mfrow=c(1,1))
lmfit <- lm(p.comp[,1]~p.comp[,2] + 0)
print(summary(lmfit))
plot(p.comp, xlab='Moment-based P-value', ylab='Permutation-based P-value', main='P-value correspondence')
lines(x = c(0,1), y = c(0,1), col='red')
legend("topleft", bty="n", legend=paste("R-squared is", format(summary(lmfit)$adj.r.squared, digits=4)))

p.comp2 <- -log10(p.comp)
par(mfrow=c(1,1))
lmfit <- lm(p.comp2[,1]~p.comp2[,2] + 0)
print(summary(lmfit))
plot(p.comp2, xlab='-log10(Moment-based P-value)', ylab='-log10(Permutation-based P-value)', main='P-value correspondence')
lines(x = c(0,10), y = c(0,10), col='red')
legend("topleft", bty="n", legend=paste("R-squared is", format(summary(lmfit)$adj.r.squared, digits=4)))
