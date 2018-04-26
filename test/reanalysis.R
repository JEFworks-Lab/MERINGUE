# Reanalysis
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

library(MUDAN)
mat <- normalizeCounts(t(as.matrix(cd)),
                       verbose=FALSE)
matnorm.info <- normalizeVariance(mat,
                                  details=TRUE,
                                  verbose=FALSE,
                                  plot=TRUE)
## log transform
matnorm <- log10(matnorm.info$mat+1)
## 30 PCs on overdispersed genes
pcs <- getPcs(matnorm[matnorm.info$ods,],
              nGenes=length(matnorm.info$ods),
              nPcs=30,
              verbose=FALSE)
## get tSNE embedding on PCs
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=10,
                    num_threads=parallel::detectCores(),
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
com <- getComMembership(pcs,
                        k=30, method=igraph::cluster_infomap,
                        verbose=FALSE)


par(mfrow=c(1,2))
col <- rainbow(length(levels(com)))[com]
names(col) <- names(com)
plotEmbedding(emb, colors=col,
              xlab=NA, ylab=NA,
              verbose=FALSE)
plot(pos$x, pos$y, col=col, pch=16, cex=2)


# merge many adjs
adjList <- lapply(3:9, function(k) {
  adj <- getAdj(pos, k=k)
})
adj <- Reduce("+", adjList) / length(adjList)
plotNetwork(pos, adj, line.power=10)

# calculate Moran's I
results <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
  value <- mat[i,]
  ape::Moran.I(value, adj)
}, mc.cores=parallel::detectCores()-1))
rownames(results) <- rownames(mat)
results <- as.data.frame(results)

results$p.adj <- stats::p.adjust(results$p.value)
results <- results[order(results$p.adj),]

#vi <- results$p.adj < 0.05
#results <- results[vi,]

g <- rownames(results)[1]
plot(pos[,1], pos[,2], col=map2col(matnorm[g,]), pch=16, cex=2)


############# Compare with SpatialDE
spatialde.results <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_final_results.csv', row.names=1)
rownames(spatialde.results) <- spatialde.results$g

a = -log10(unlist(results$p.adj))
a[is.infinite(a)] <- NA
b = -log10(spatialde.results[rownames(results),]$qval)
b[is.infinite(b)] <- NA
plot(a, b,
     ylab = "-log10(p-value) for SpatialDE",
     xlab = "-log10(p-value) for MERingue")
cor.test(a,b, use="pairwise.complete.obs")

vi <- a>-log10(0.05) & b>-log10(0.05)
plot(a[vi], b[vi],
     ylab = "-log10(p-value) for SpatialDE",
     xlab = "-log10(p-value) for MERingue")


############# Compare with trensceeek
library('trendsceek')

pp = pos2pp(pos[,1:2])
log.fcn = log10
pp = set_marks(pp, t(matnorm['Ptn',]), log.fcn = log.fcn)

pp2plot = pp_select(pp)

##set parameters
nrand = 100
ncores = 1
##run
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

trendstat_list[['supstats_wide']]

g <- 'Ptn'
plot(pos[,1], pos[,2], col=map2col(matnorm[g,]), pch=16, cex=3)
