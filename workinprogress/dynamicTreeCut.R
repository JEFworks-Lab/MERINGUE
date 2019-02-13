source('../R/main.R')
source('../R/process.R')
source('../R/helper.R')

#install.packages('dynamicTreeCut')
library(dynamicTreeCut)

# Determine gene clusters using dynamicTreeCut

# mouse olfactory bulb
cd <- read.csv('../../../SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('../../../SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# Clean
counts <- cleanCounts(t(cd), min.reads=10, min.detected=10)
mat <- normalizeCounts(counts)

# Analyze
adj <- getSpatialWeights(pos, plot=FALSE)
results <- getSpatialPatterns(mat, adj)
# Get significant hits
vi <- results$p.adj < 0.01
results.sig <- rownames(results)[vi]
results.sig

# Cluster on sig genes
m <- as.matrix(mat[results.sig,])
m <- apply(m, 1, winsorize) # get rid of outliers
m <- scale(m)
d <- dist(t(m))
hc <- hclust(d)
par(mfrow=c(1,1))
plot(hc)
# static tree cut
table(cutree(hc, k=5))
# dynamic tree cut
groups <- cutreeDynamic(dendro=hc, distM=as.matrix(d), method='hybrid', minClusterSize=10, deepSplit=0)
names(groups) <- hc$labels
groups <- factor(groups)
table(groups)

par(mfrow=c(length(levels(groups)), 2), mar=rep(1,4))
prs <- lapply(levels(groups), function(g) {
    # summarize as first pc if more than 1 gene in group
    if(sum(groups==g)>1) {
        m <- winsorize(mat[results.sig[groups==g],])
        ps <- colMeans(m)
        pc <- prcomp(m)
        pr <- pc$rotation[,1]
        # double check direction is same
        if(cor(ps, pr)<0) {
            pr <- -pr
        }
    } else {
        pr <- winsorize(mat[results.sig[groups==g],])
    }
    interpolate(pos, pr, main=paste0("Pattern ", g, " : ", sum(groups==g), " genes"), plot=TRUE)
    return(pr)
})
names(prs) <- levels(groups)

gsub <- names(groups)[groups==levels(groups)[6]][1:10]
par(mfrow=c(length(gsub), 2), mar=rep(1,4))
prs <- lapply(gsub, function(g) {
    pr <- winsorize(mat[g,])
    interpolate(pos, pr, main=g, plot=TRUE)
    return(pr)
})
