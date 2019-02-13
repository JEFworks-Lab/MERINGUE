# Parse data
parseData <- function(jsfile) {
  library(jsonlite)
  result <- jsonlite::fromJSON(jsfile, flatten=TRUE)

  pos <- data.frame(result$barcode, x=result$x, y=result$y)
  vi <- !duplicated(pos)
  table(vi)
  pos <- pos[vi,]
  rownames(pos) <- pos[,1]
  pos <- pos[,-1]

  dat <- data.frame(cell=result$barcode, gene=result$gene, value=result$hits)
  library(reshape2)
  cd <- dcast(data = dat, formula = cell~gene, fun.aggregate = sum, value.var = "value")
  rownames(cd) <- cd[,1]
  cd <- cd[,-1]
  cd <- t(cd)

  return(list(cd=cd, pos=pos))
}
BCL1 <- parseData('data/BC_layer1.json_.gz')
BCL2 <- parseData('data/BC_layer2.json_.gz')
BCL3 <- parseData('data/BC_layer3.json_.gz')
BCL4 <- parseData('data/BC_layer4.json_.gz')


# Align layers
g <- 'ACTB'
par(mfrow=c(4,2))
p1 <- BCL1$pos
p1[,1] <- p1[,1] - mean(p1[,1])
p1[,2] <- p1[,2] - mean(p1[,2])
interpolate(p1, BCL1$cd[g,], main='Layer 1')
p2 <- rotatePos(BCL2$pos, (-12*pi)/180)
p2[,1] <- p2[,1] - mean(p2[,1])
p2[,2] <- p2[,2] - mean(p2[,2])
interpolate(p2, BCL2$cd[g,], main='Layer 2')
p3 <- rotatePos(BCL3$pos, (-148*pi)/180)
p3[,1] <- p3[,1] - mean(p3[,1])
p3[,2] <- p3[,2] - mean(p3[,2])
interpolate(p3, BCL3$cd[g,], main='Layer 3')
p4 <- rotatePos(BCL4$pos, (-175*pi)/180)
p4[,1] <- p4[,1] - mean(p4[,1])
p4[,2] <- p4[,2] - mean(p4[,2])
interpolate(p4, BCL4$cd[g,], main='Layer 4')
plot(p4)
points(p3)
points(p2)
points(p1)

BCL1$pos <- p1
BCL2$pos <- p2
BCL3$pos <- p3
BCL4$pos <- p4


################################# transcriptional clustering

genes.have <- Reduce(intersect, list(
  rownames(BCL1$cd),
  rownames(BCL2$cd),
  rownames(BCL3$cd),
  rownames(BCL4$cd)

))
counts <- cbind(
  BCL1$cd[genes.have,],
  BCL2$cd[genes.have,],
  BCL3$cd[genes.have,],
  BCL4$cd[genes.have,]
)
colnames(counts) <- c(
  paste0('L1-', colnames(BCL1$cd)),
  paste0('L2-', colnames(BCL2$cd)),
  paste0('L3-', colnames(BCL3$cd)),
  paste0('L4-', colnames(BCL4$cd))
)
layer <-  c(
  rep('L1', ncol(BCL1$cd)),
  rep('L2', ncol(BCL2$cd)),
  rep('L3', ncol(BCL3$cd)),
  rep('L4', ncol(BCL4$cd))
)
names(layer) <- colnames(counts)

cc <- cleanCounts(counts, min.detected=100)
hist(log10(colSums(as.matrix(cc))+1))
pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
plotEmbedding(pos1, colors=colSums(as.matrix(cc)), cex=2, main='Lib size')
plotEmbedding(pos1, colors=colSums(as.matrix(cc)>0), cex=2, main='Lib complexity')

mat <- normalizeCounts(cc)
plotEmbedding(pos1, colors=colSums(as.matrix(mat)), cex=2, main='Lib size')
plotEmbedding(pos1, colors=colSums(as.matrix(mat)>0), cex=2, main='Lib complexity')

library(MUDAN)
ods <- normalizeVariance(cc, plot=TRUE, details=TRUE, alpha=0.01)
#m <- ods$mat[intersect(sig.genes, rownames(ods$mat)[ods$ods]),]
#m <- ods$mat[ods$ods,]
m <- mat[ods$ods,]
dim(m)
pcs <- getPcs(m, nPcs=100, nGenes=nrow(m))
d <- dist(pcs, method='euc')
emb <- Rtsne::Rtsne(d,
                    is_distance=TRUE,
                    perplexity=30,
                    num_threads=2,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)

par(mfrow=c(1,2), mar=rep(2,4))
plotEmbedding(emb, groups=layer, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='topright', main='Layers', cex=1)

com <- getComMembership(pcs,
                        k=30, method=igraph::cluster_walktrap,
                        verbose=FALSE)
plotEmbedding(emb, groups=com, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='bottomleft')


###################### look at differential genes
dg <- getDifferentialGenes(cc, com)
dg.sig <- lapply(dg, function(x) {
  x$pv.adj <- p.adjust(pnorm(-x$Z), n=nrow(counts)*ncol(counts))
  x <- na.omit(x)
  #x <- x[x$highest,]
  #x <- x[abs(x$Z)>1.96,]
  x <- x[x$Z>1.96,]
  x[order(x$Z, decreasing=TRUE),]
})
names(dg.sig) <- names(dg)
lapply(dg.sig, nrow)
dg.sig[[3]]

pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
plotEmbedding(pos1, groups=com, cex=2)
g <- 'PPARG'
gexp <- mat[g,]
hist(gexp)
plotEmbedding(pos1, colors=gexp, cex=2, main=g)


##################### Annotate

annot <- as.character(com)
names(annot) <- names(com)
annot[com==1] <- 'INV'
annot[com==2] <- 'Stroma'
annot[com==3] <- 'Adipose'
annot[com==4] <- 'DCIS'
plotEmbedding(emb, groups=annot, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='topright', main='Clusters', cex=1)

# Visualize
par(mfrow=c(2,2))
library(jpeg)

pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
ima <- readJPEG("data/HE_layer1_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1, main='Layer 1')
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[rownames(pos1),], col=fac2col(annot[rownames(pos1)], v=0.8), cex=1, lwd=3)

pos1 <- BCL2$pos
rownames(pos1) <- paste0('L2-', rownames(pos1))
ima <- readJPEG("data/HE_layer2_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1, main='Layer 2')
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[rownames(pos1),], col=fac2col(annot[rownames(pos1)], v=0.8), cex=1, lwd=3)

pos1 <- BCL3$pos
rownames(pos1) <- paste0('L3-', rownames(pos1))
ima <- readJPEG("data/HE_layer3_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1, main='Layer 3')
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[rownames(pos1),], col=fac2col(annot[rownames(pos1)], v=0.8), cex=1, lwd=3)

pos1 <- BCL4$pos
rownames(pos1) <- paste0('L4-', rownames(pos1))
ima <- readJPEG("data/HE_layer4_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1, main='Layer 4')
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[rownames(pos1),], col=fac2col(annot[rownames(pos1)], v=0.8), cex=1, lwd=3)


# Restrict to INV and ductal breast cancer cells
#subset <- names(annot)[!(annot %in% c('Stroma', 'Epithelium', 'Immune Infiltrate'))]
#subset <- names(annot)[annot %in% c('DCIS', 'INV')]
subset <- names(com)[com %in% c(1,2,4)]

par(mfrow=c(2,2))
library(jpeg)

pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
ima <- readJPEG("data/HE_layer1_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)

pos1 <- BCL2$pos
rownames(pos1) <- paste0('L2-', rownames(pos1))
ima <- readJPEG("data/HE_layer2_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)

pos1 <- BCL3$pos
rownames(pos1) <- paste0('L3-', rownames(pos1))
ima <- readJPEG("data/HE_layer3_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)

pos1 <- BCL4$pos
rownames(pos1) <- paste0('L4-', rownames(pos1))
ima <- readJPEG("data/HE_layer4_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)



##################### Get highly spatially variable genes (union)
getI <- function(BCL1, subset) {
  w <- getSpatialWeights(BCL1$pos[subset,], klist=3)
  plotNetwork(pos=BCL1$pos[subset,], adj=w)
  mat <- normalizeCounts(cleanCounts(BCL1$cd[,subset], min.reads = 10))
  I <- getSpatialPatterns(mat, w)
  results <- I
  vi <- results$p.adj < 0.05
  #vi <- results$p.value < 0.05
  vi[is.na(vi)] <- FALSE
  results.sig <- rownames(results)[vi]
  return(results.sig)
}

par(mfrow=c(2,2))
subset1 <- gsub('L1-', '', subset[grepl('L1', subset)])
sig1 <- getI(BCL1, subset1)
subset2 <- gsub('L2-', '', subset[grepl('L2', subset)])
sig2 <- getI(BCL2, subset2)
subset3 <- gsub('L3-', '', subset[grepl('L3', subset)])
sig3 <- getI(BCL3, subset3)
subset4 <- gsub('L4-', '', subset[grepl('L4', subset)])
sig4 <- getI(BCL4, subset4)





##### Use 3D information
pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
pos1 <- pos1[intersect(subset, rownames(pos1)),]

pos2 <- BCL2$pos
rownames(pos2) <- paste0('L2-', rownames(pos2))
pos2 <- pos2[intersect(subset, rownames(pos2)),]

pos3 <- BCL3$pos
rownames(pos3) <- paste0('L3-', rownames(pos3))
pos3 <- pos3[intersect(subset, rownames(pos3)),]

pos4 <- BCL4$pos
rownames(pos4) <- paste0('L4-', rownames(pos4))
pos4 <- pos4[intersect(subset, rownames(pos4)),]

x <- c(pos1[,1], pos2[,1], pos3[,1], pos4[,1])
y <- c(pos1[,2], pos2[,2], pos3[,2], pos4[,2])
z <- c(rep(1, nrow(pos1)),
       rep(2, nrow(pos2)),
       rep(3, nrow(pos3)),
       rep(4, nrow(pos4))
)
group <- c(rep('L1', nrow(pos1)),
           rep('L2', nrow(pos2)),
           rep('L3', nrow(pos3)),
           rep('L4', nrow(pos4))
)
names(group) <- names(x) <- names(y) <- names(z) <- c(rownames(pos1), rownames(pos2), rownames(pos3), rownames(pos4))
library(rgl)
plot3d(x,y,z, col=rainbow(4)[as.factor(group)], size=10)

ctA <- rownames(pos1)
ctB <- rownames(pos2)
pos <- rbind(pos1, pos2)
pos <- as.matrix(pos)
w1 <- getMnn(ctA, ctB, pos=pos, k=3)

ctA <- rownames(pos2)
ctB <- rownames(pos3)
pos <- rbind(pos2, pos3)
pos <- as.matrix(pos)
w2 <- getMnn(ctA, ctB, pos=pos, k=3)

ctA <- rownames(pos3)
ctB <- rownames(pos4)
pos <- rbind(pos3, pos4)
pos <- as.matrix(pos)
w3 <- getMnn(ctA, ctB, pos=pos, k=3)

w <- matrix(0, nrow=length(group), ncol=length(group))
dim(w)
rownames(w) <- colnames(w) <- names(group)
w[rownames(w1), colnames(w1)] <- w1
w[rownames(w2), colnames(w2)] <- w2
w[rownames(w3), colnames(w3)] <- w3

plot3d(x[names(group)],y[names(group)],z[names(group)], col=fac2col(group), size=10)
idx <- which(w>0, arr.ind = T)
for(i in seq_len(nrow(idx))) {
  n1 <- rownames(w)[idx[i,1]]
  n2 <- rownames(w)[idx[i,2]]
  segments3d(
    x=c(x[n1], x[n2]),
    y=c(y[n1], y[n2]),
    z=c(z[n1], z[n2]),
    col='grey'
  )
}

mat <- normalizeCounts(cleanCounts(counts[, subset], min.reads = 10))
I <- getSpatialPatterns(mat, w)
results <- I
vi <- results$p.adj < 0.2
#vi <- results$p.value < 0.05
vi[is.na(vi)] <- FALSE
results.sig <- rownames(results)[vi]
results.sig


library(UpSetR)
sig <- list('L1'=sig1, 'L2'=sig2, 'L3'=sig3, 'L4'=sig4, 'all'=results.sig)
upset(fromList(sig))

sort(Reduce(intersect, sig))
intersect(sig1, sig2)
intersect(sig2, sig3)
intersect(sig3, sig4)

g =  "FN1"
gexp <- scale(mat[g,intersect(subset, names(group))])[,1]
library(rgl)
plot3d(x[subset],y[subset],z[subset], col=map2col(winsorize(gexp[subset]), colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10, main=g)
planes3d(0, 0, 1, -1, alpha = 0.1)
planes3d(0, 0, 1, -2, alpha = 0.1)
planes3d(0, 0, 1, -3, alpha = 0.1)
planes3d(0, 0, 1, -4, alpha = 0.1)

## plot expression of gene
par(mfrow=c(2,2))
subset1 <- gsub('L1-', '', subset[grepl('L1', subset)])
pos1 <- BCL1$pos
ima <- readJPEG("data/HE_layer1_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(winsorize(mat[g,subset[grepl('L1', subset)]])), cex=2, lwd=3, pch=16)

subset1 <- gsub('L2-', '', subset[grepl('L2', subset)])
pos1 <- BCL2$pos
ima <- readJPEG("data/HE_layer2_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(winsorize(mat[g,subset[grepl('L2', subset)]])), cex=2, lwd=3, pch=16)

subset1 <- gsub('L3-', '', subset[grepl('L3', subset)])
pos1 <- BCL3$pos
ima <- readJPEG("data/HE_layer3_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(winsorize(mat[g,subset[grepl('L3', subset)]])), cex=2, lwd=3, pch=16)

subset1 <- gsub('L4-', '', subset[grepl('L4', subset)])
pos1 <- BCL4$pos
ima <- readJPEG("data/HE_layer4_BC.jpg")
ima <- ima[,,1] # greyscale via red channel
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(winsorize(mat[g,subset[grepl('L4', subset)]])), cex=2, lwd=3, pch=16)

