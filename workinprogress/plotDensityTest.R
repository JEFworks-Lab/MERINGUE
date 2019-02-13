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

# Get transcriptional clusters
library(MUDAN)
matnorm <- normalizeVariance(counts, plot=TRUE, details = TRUE, alpha=0.2)
pcs <- getPcs(mat[matnorm$ods,], nPcs=10)

d <- dist(pcs, method='euc')
emb <- Rtsne::Rtsne(d,
                    is_distance=TRUE,
                    perplexity=10,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)

com <- getComMembership(pcs,
                        k=30, method=igraph::cluster_louvain,
                        verbose=FALSE)

annot <- as.character(com); names(annot) <- names(com)
annot[com==1] <- '3: Glomerular Layer'
annot[com==2] <- '1: Granular Cell Layer'
annot[com==3] <- '2: Mitral Cell Layer'
annot[com==4] <- '4: Olfactory Nerve Layer'
annot <- as.factor(annot)

par(mfrow=c(1,2), mar=rep(2,4))
plotEmbedding(emb, groups=annot, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='bottomleft')
plotEmbedding(pos, groups=annot, cex=1, xlab=NA, ylab=NA)

##### Test density plots
par(mfrow=c(1,2))
plotDensity(pos)
plotDensitySubset(pos, names(com)[com==1])

