# get neighbor's neighbor
set.seed(0)
library(MERingue)
data(mOB)

# Clean and normalize data
pos <- mOB$pos
cd <- mOB$counts

counts <- cleanCounts(cd, min.reads=10, min.lib.size=10)
pos <- pos[colnames(counts),]
mat <- normalizeCounts(counts, log=FALSE)

w <- voronoiAdjacency(pos, filterDist = 2)
par(mfrow=c(1,1), mar=rep(5,4))
plotNetwork(pos, w)

w0 <- w
invisible(lapply(seq_len(nrow(w)), function(i) {
  nn <- colnames(w)[which(w[which(w[i,]==1),]==1, arr.ind=TRUE)[,2]]
  w0[i,nn] <<- 1/2
}))

table(w0==w)
plotNetwork(pos, w0)

I <- getSpatialPatterns(mat, w0)
results.filter.df <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = w0,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.1,
                                        verbose = TRUE,
                                        details=TRUE)

results.filter <- rownames(results.filter.df)
length(results.filter)

g <- results.filter[2]
gexp <- mat[g,]
#gexp <- log10(mat[g,]+1)
moranTest(gexp, w)
hist(gexp)
plot(pos, col=MERingue:::map2col(gexp), pch=16)

## spatial cross correlation
scc <- spatialCrossCorMatrix(as.matrix(mat[results.filter,]), w0)
ggroup <- groupSigSpatialPatterns(pos, as.matrix(mat[results.filter,]), scc, power=1, hclustMethod='ward.D', deepSplit = 3)
gcol <- rainbow(length(levels(ggroup$groups)), v=0.5)[ggroup$groups]
names(gcol) <- names(ggroup$groups)

## compare
diffgexp <- dg.sig
#diffgexp <- lapply(levels(anova.groups), function(x) names(anova.groups)[anova.groups==x])
#names(diffgexp) <- levels(anova.groups)
spatgexp <- lapply(levels(ggroup$groups), function(x) {
  names(ggroup$groups)[ggroup$groups==x]
})
names(spatgexp) <- paste0('spatial', levels(ggroup$groups))

library(UpSetR)
upset(UpSetR::fromList(c(diffgexp, spatgexp)), sets=names(c(diffgexp, spatgexp)), keep.order=TRUE, order.by="degree")

setdiff(unlist(spatgexp[[6]]), unlist(diffgexp))
g <- 'Apoe'
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))
hist(log10(mat[g,]+1))


## may just be markign multiple subpops, check anova
pv <- sapply(1:nrow(mat), function(i) {
  mydataframe <- data.frame(y=mat[i,], ig=annot)
  fit <- aov(y ~ ig, data=mydataframe)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
names(pv) <- rownames(mat)
pv.sig <- names(pv)[p.adjust(pv) < 0.05] ## bonferonni
anova <- pv.sig

setdiff(unlist(spatgexp[[1]]), anova)
g <- 'Apoe'
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))
hist(log10(mat[g,]+1))
