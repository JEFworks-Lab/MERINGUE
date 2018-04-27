# Reanalysis 3
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# Clean
counts <- cleanCounts(t(cd), min.reads=10, min.detected=10)
mat <- normalizeCounts(counts)

par(mfrow=c(1,1))
adj <- getSpatialWeights(pos, plot=TRUE)
results <- getSpatialPatterns(mat, adj)
head(results)

vi <- results$p.adj < 0.05
table(vi)
results.sig <- rownames(results)[vi]
head(results.sig)

g <- results.sig[6]
interpolate(pos, mat[g,])

cv <- cor(t(as.matrix(mat[results.sig,])))
hc <- hclust(as.dist(1-cv))
heatmap(cv[hc$labels, hc$labels], Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorRampPalette(c('blue', 'white', 'red'))(100), labRow=NA, labCol=NA)
groups <- cutree(hc, 5)

pc <- prcomp(mat[results.sig[groups==4],])
pr <- pc$rotation[,1]
interpolate(pos, pr)
