## Compare different ways to characterize gene expression correlations to identify patterns

# Sample analysis using OB data
source('R/nearestNeighbor.R')
source('R/correlationStatistics.R')
source('R/main.R')
source('R/process.R')
source('R/plot.R')

# Data
cd <- read.csv("~/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv", row.names=1)
pos <- read.csv("~/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv", row.names=1)
cd <- as.matrix(cd[rownames(pos),])
sample_info <- pos
pos <- pos[,1:2]
devtools::use_data(cd, pos)

# Process
counts <- cleanCounts(t(as.matrix(cd)), min.reads=100, min.detected=100)
mat <- normalizeCounts(counts)

# Regular clustering
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

# Analyze spatial correlations
w <- getSpatialWeights(pos, klist=3)
par(mfrow=c(1,1))
plotNetwork(pos, w)
I <- getSpatialPatterns(mat, w, verbose=TRUE, ncores=2)
results <- I
vi <- results$p.adj < 0.05
vi[is.na(vi)] <- FALSE
table(vi)
results.sig <- rownames(results)[vi]

## regular correlation
rc <- cor(t(as.matrix(mat[results.sig,])))
d <- as.dist(1-rc)
rc.results <- groupSigSpatialPatterns(pos, mat, d, minClusterSize=10)

library(gplots)
library(RColorBrewer)
m <- mat[names(sort(rc.results$groups)),names(sort(annot))]
m <- t(scale(t(m)))
m[m < -2] <- -2
m[m > 2] <- 2
heatmap.2(trace="none", m, Rowv=NA, Colv=NA,
          ColSideColors = rainbow(length(levels(annot)))[annot[colnames(m)]],
          RowSideColors = rainbow(length(levels(rc.results$groups)))[rc.results$groups[rownames(m)]],
          scale="row", col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))


## Spatial cross correlation
scc <- do.call(cbind, parallel::mclapply(results.sig, function(g1) {
  print(g1)
  unlist(lapply(results.sig, function(g2) {
    x <- mat[g1,]
    y <- mat[g2,]
    spatialIntraCrossCor(x,y,w)
  }))
}, mc.cores=2))
dim(scc)
rownames(scc) <- colnames(scc) <- results.sig
d <- as.dist((-scc - min(-scc)))
scc.results <- groupSigSpatialPatterns(pos, mat, d, minClusterSize=10)

m <- mat[names(sort(scc.results$groups)),names(sort(annot))]
m <- t(scale(t(m)))
m[m < -2] <- -2
m[m > 2] <- 2
heatmap.2(trace="none", m, Rowv=NA, Colv=NA,
          ColSideColors = rainbow(length(levels(annot)))[annot[colnames(m)]],
          RowSideColors = rainbow(length(levels(scc.results$groups)))[scc.results$groups[rownames(m)]],
          scale="row", col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))


## Identify pair of genes with low correlation but high spatial cross correlation
table(rc.results$groups, scc.results$groups)
par(mfrow=c(2,1), mar=rep(5,4))
plot(as.dendrogram(rc.results$hc)[[1]])
plot(as.dendrogram(scc.results$hc)[[1]])

sort(sapply(results.sig, function(i) lm(scc[i,]~rc[i,])$coefficients[2]), decreasing=TRUE)

i <- 'Cdhr1'
plot(rc[i,], scc[i,])
text(rc[i,], scc[i,], labels=results.sig)
j <- 'Gabra1'

rc[i,j]
scc[i,j]

par(mfrow=c(2,2))
interpolate(pos, mat[i,], main=i, zlim=c(-2,2))
interpolate(pos, mat[j,], main=j, zlim=c(-2,2))

## weighted correlation by LISA
lisa.pv <- do.call(cbind, lapply(results.sig, function(g) {
  Ii <- lisaTest(mat[g,], w)
  pv <- -log10(Ii$p.value)
  names(pv) <- rownames(Ii)
  return(pv)
}))
colnames(lisa.pv) <- results.sig



lisac <- cor(lisa.pv)
d <- as.dist(1-lisac)
scc.results <- groupSigSpatialPatterns(pos, mat, d)
