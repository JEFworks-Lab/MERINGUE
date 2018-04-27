# Within cell type heterogeneity

# approximate using
# mouse olfactory bulb
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# Clean
counts <- cleanCounts(t(cd), min.reads=10, min.detected=10)
mat <- normalizeCounts(counts)

# Analyze
adj <- getSpatialWeights(pos, plot=TRUE)
results <- getSpatialPatterns(mat, adj)
results.sig <- rownames(results)[results$p.adj < 0.05]

# Cluster cells
library(MUDAN)
pcs <- getPcs(mat[results.sig,],
              nGenes=ncol(mat),
              nPcs=30,
              verbose=FALSE)
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=5,
                    num_threads=parallel::detectCores(),
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
com <- getComMembership(pcs,
                        k=50, method=igraph::cluster_infomap,
                        verbose=FALSE)
par(mfrow=c(1,2))
plotEmbedding(emb, com, xlab=NA, ylab=NA,
              mark.clusters=TRUE, alpha=0.5, mark.cluster.cex=0.5,
              verbose=FALSE) ## plot
plot(pos, col=map2col(as.numeric(com)), pch=16)

# subset
mat.sub <- mat[, com==2]
pos.sub <- pos[com==2,]
# Analyze
adj.sub <- getSpatialWeights(pos.sub, plot=TRUE)
results.sub <- getSpatialPatterns(mat.sub, adj.sub)

groups <- groupSigSpatialPatterns(pos.sub, mat.sub, results.sub, k=5)
par(mfrow=c(length(levels(groups)), 2), mar=rep(1,4))
prs <- lapply(levels(groups), function(g) {
  # summarize as first pc if more than 1 gene in group
  if(sum(groups==g)>1) {
    pc <- prcomp(mat[results.sig[groups==g],])
    pr <- pc$rotation[,1]
  } else {
    pr <- mat[results.sig[groups==g],]
  }
  interpolate(pos, pr, main=paste0("Pattern ", g, " : ", sum(groups==g), " genes"), plot=TRUE)
  return(pr)
})
names(prs) <- levels(groups)

