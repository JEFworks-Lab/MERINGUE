# Reanalysis 2
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

par(mfrow=c(1,1))
adj <- getSpatialWeights(pos, plot=TRUE)

# test sensitivity to different normalizations
results1 <- getSpatialPatterns(log10(t(cd)+1), adj)
head(results1)

# check that genes not expressed are generally not significant
gexp <- log10(colSums(cd>0)+1)
hist(gexp)
test <- names(which(gexp < 1))
hist(results1[test,]$p.adj)

# CPM normalization
counts <- t(cd)
mat <- Matrix::t(Matrix::t(counts)/Matrix::colSums(counts))
mat <- mat*1e6
mat <- log10(mat+1)
results2 <- getSpatialPatterns(mat, adj)

head(results1)
head(results2)

a <- -log10(results1$p.adj)
b <- -log10(results2$p.adj)
names(a) <- rownames(results1)
names(b) <- rownames(results2)
a[is.infinite(a)] <- NA
b[is.infinite(b)] <- NA
vi <- b>1 | a>1
table(vi)
plot(a[vi], b[vi], pch=16)

g <- rownames(results1)[4]
g <- rownames(results2)[11]
plot(pos, col=map2col(mat[g,]), pch=16)

# Conclusion, definitely need proper normalization to ensure autocorrelations are appropriate
