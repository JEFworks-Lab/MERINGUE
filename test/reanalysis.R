# Reanalysis

cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]

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
plotEmbedding(emb, com,
              xlab=NA, ylab=NA,
              mark.clusters=TRUE, alpha=0.5, mark.cluster.cex=0.5,
              verbose=FALSE)


par(mfrow=c(1,1))
plot(pos[,1], pos[,2], col=com, pch=16, cex=3)




# plot
adj <- getAdj(pos[,1:2], k=3)
plotNetwork(pos[,1:2], adj)

# calculate Moran's I
results <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
  value <- mat[i,]
  ape::Moran.I(value, adj)
}, mc.cores=parallel::detectCores()-1))
rownames(results) <- rownames(mat)
results <- as.data.frame(results)

results$p.adj <- stats::p.adjust(results$p.value)
results <- results[order(results$p.adj),]

head(results)

g <- rownames(results)[1]
plot(pos[,1], pos[,2], col=map2col(matnorm[g,]), pch=16, cex=3)




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
