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
p2 <- rotatePos(BCL2$pos, (-10*pi)/180)
p2[,1] <- p2[,1] - mean(p2[,1])
p2[,2] <- p2[,2] - mean(p2[,2])
interpolate(p2, BCL2$cd[g,], main='Layer 2')
p3 <- rotatePos(BCL3$pos, (-150*pi)/180)
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

library(MUDAN)
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

mat <- normalizeCounts(counts)
ods <- normalizeVariance(counts, plot=TRUE, details=TRUE)
#m <- ods$mat[intersect(sig.genes, rownames(ods$mat)),]
m <- ods$mat[ods$ods,]
dim(m)
pcs <- getPcs(m, nPcs=30, nGenes=nrow(m))
d <- dist(pcs, method='euc')
emb <- Rtsne::Rtsne(d,
                    is_distance=TRUE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)

par(mfrow=c(1,1), mar=rep(2,4))
plotEmbedding(emb, groups=layer, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='bottomleft')

com <- getComMembership(pcs,
                        k=30, method=igraph::cluster_infomap,
                        verbose=FALSE)

annot <- as.character(com)
names(annot) <- names(com)
annot[com==1] <- 'Stroma'
annot[com==2] <- 'Epithelium'
annot[com==3] <- 'Ductal Carcinoma'
annot[com==4] <- 'Ductal Carcinoma + Stroma'
annot[com==7] <- 'INV'
annot[com==5] <- 'INV + Stroma'
annot[com==6] <- 'Immune Infiltrate'
plotEmbedding(emb, groups=annot, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='bottomleft')


subset <- names(annot)[annot %in% c('Ductal Carcinoma', 'INV')]


library(jpeg)
pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
ima <- readJPEG("data/HE_layer1_BC.jpg")
plotEmbedding(pos1, groups=annot[rownames(pos1)], cex=2, show.legend=TRUE, legend.x='bottomleft')
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)


library(jpeg)
pos1 <- BCL2$pos
rownames(pos1) <- paste0('L2-', rownames(pos1))
ima <- readJPEG("data/HE_layer2_BC.jpg")
plotEmbedding(pos1, groups=annot[rownames(pos1)], cex=2, show.legend=TRUE, legend.x='bottomleft')
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)


library(jpeg)
pos1 <- BCL3$pos
rownames(pos1) <- paste0('L3-', rownames(pos1))
ima <- readJPEG("data/HE_layer3_BC.jpg")
plotEmbedding(pos1, groups=annot[rownames(pos1)], cex=2, show.legend=TRUE, legend.x='bottomleft')
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)


library(jpeg)
pos1 <- BCL4$pos
rownames(pos1) <- paste0('L4-', rownames(pos1))
ima <- readJPEG("data/HE_layer4_BC.jpg")
plotEmbedding(pos1, groups=annot[rownames(pos1)], cex=2, show.legend=TRUE, legend.x='bottomleft')
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(subset, rownames(pos1)),], col=fac2col(annot[intersect(subset, rownames(pos1))], v=0.8), cex=2, lwd=3)


################################ look at differential genes
dg <- getDifferentialGenes(mat, annot)
dg.sig <- lapply(dg, function(x) {
  x$pv.adj <- p.adjust(pnorm(-x$Z), n=nrow(counts)*ncol(counts))
  x <- na.omit(x)
  x <- x[x$highest,]
  #x <- x[abs(x$Z)>1.96,]
  x <- x[x$Z>1.96,]
  x[order(x$Z, decreasing=TRUE),]
})
names(dg.sig) <- names(dg)
lapply(dg.sig, nrow)
dg.sig[['Cancer-Ductal']]
dg.sig[['Cancer-INV']]

pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
g <- 'TNC'
gexp <- counts[g,]
plotEmbedding(pos1, colors=gexp, cex=2, main=g)




# Get highly spatially variable genes (union)
getI <- function(BCL1, subset) {
  w <- getSpatialWeights(BCL1$pos[subset,], klist=3)
  plotNetwork(pos=BCL1$pos[subset,], adj=w)
  mat <- normalizeCounts(BCL1$cd[,subset])
  I <- getSpatialPatterns(mat, w)
  results <- I
  vi <- results$p.adj < 0.05
  #vi <- results$p.value < 0.05
  vi[is.na(vi)] <- FALSE
  results.sig <- rownames(results)[vi]
  return(results.sig)
}

subset1 <- gsub('L1-', '', subset[grepl('L1', subset)])
sig1 <- getI(BCL1, subset1)
sig1

par(mfrow=c(2,2))
g <- "TAGLN2"
g <- "CAMP"
g <- "SAA1"
interpolate(BCL1$pos, BCL1$cd[g,], main=g)
interpolate(BCL1$pos, BCL1$cd[g,subset1], main=g)

subset1 <- gsub('L1-', '', subset[grepl('L1', subset)])
library(jpeg)
pos1 <- BCL1$pos
ima <- readJPEG("data/HE_layer1_BC.jpg")
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(log10(BCL1$cd[g,subset1]+1)), cex=2, lwd=3, pch=16)

subset1 <- gsub('L2-', '', subset[grepl('L2', subset)])
library(jpeg)
pos1 <- BCL2$pos
ima <- readJPEG("data/HE_layer2_BC.jpg")
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(log10(BCL2$cd[g,subset1]+1)), cex=2, lwd=3, pch=16)

subset1 <- gsub('L3-', '', subset[grepl('L3', subset)])
library(jpeg)
pos1 <- BCL3$pos
ima <- readJPEG("data/HE_layer3_BC.jpg")
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(log10(BCL3$cd[g,subset1]+1)), cex=2, lwd=3, pch=16)

subset1 <- gsub('L4-', '', subset[grepl('L4', subset)])
library(jpeg)
pos1 <- BCL4$pos
ima <- readJPEG("data/HE_layer4_BC.jpg")
par(mfrow=c(1,1))
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[subset1,], col=map2col(log10(BCL4$cd[g,subset1]+1)), cex=2, lwd=3, pch=16)










sig2 <- getI(BCL2, subset)
sig3 <- getI(BCL3)
sig4 <- getI(BCL4)
sig.genes <- Reduce(union, list(sig1, sig2, sig3, sig4))
length(sig.genes)

# Get patterns
#scc <- spatialCrossCorMatrix(as.matrix(BCL1$cd[sig1,]), getSpatialWeights(BCL1$pos, klist=3))
#dim(scc)
#d <- as.dist((-scc - min(-scc))^(1/3))
#d <- as.dist(1-cor(t(as.matrix(BCL1$cd[sig1,]))))
#ggroup <- groupSigSpatialPatterns(pos=BCL1$pos, as.matrix(BCL1$cd[sig1,]), d, minClusterSize = 30)





# align by nearest neighbor
pos1 <- BCL1$pos
pos2 <- BCL4$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
rownames(pos2) <- paste0('L4-', rownames(pos2))
gexp1 <- BCL1$cd
gexp2 <- BCL4$cd
colnames(gexp1) <- paste0('L1-', colnames(gexp1))
colnames(gexp2) <- paste0('L4-', colnames(gexp2))
sig.genes <- intersect(intersect(c(sig1,sig2), rownames(BCL1$cd)), rownames(BCL4$cd))
sig.genes
gexp1 <- gexp1[sig.genes,]
gexp2 <- gexp2[sig.genes,]

ctA <- rownames(pos1)
ctB <- rownames(pos2)
pos <- rbind(pos1, pos2)
pos <- as.matrix(pos)
group <- as.factor((c(ctA, ctB) %in% ctA)*1+1)
names(group) <- c(ctA, ctB)
table(group)
gexpdf <- data.frame(t(cbind(gexp1, gexp2)))
w <- getMnn(ctA, ctB, pos=gexpdf, k=1)
plotNetwork(pos, w, col=rainbow(2)[group], line.col = 'grey')

library(rgl)
x <- c(pos1[,1], pos2[,1])
y <- c(pos1[,2], pos2[,2])
z <- c(rep(1, nrow(pos1)), rep(2, nrow(pos2)))
names(x) <- names(y) <- names(z) <- c(rownames(pos1), rownames(pos2))
plot3d(x,y,z, col=rainbow(2)[group], size=10)
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

# manually rotate
pos2 <- rotatePos(pos2, pi)
pos <- rbind(pos1, pos2)
pos <- as.matrix(pos)

# look at positional differences
lmfitx <- plotNeighborBoxplot(pos[,1], pos[,1], ctA, ctB, w)
dx <- lmfitx$coefficients[1]
dx
lmfity <- plotNeighborBoxplot(pos[,2], pos[,2], ctA, ctB, w)
dy <- lmfity$coefficients[1]
dy
pos2[,1] <- pos2[,1]-dx
pos2[,2] <- pos2[,2]-dy

lmfitx <- plotNeighborBoxplot(pos[,1], pos[,1], ctB, ctA, w)
dx <- lmfitx$coefficients[1]
dx
lmfity <- plotNeighborBoxplot(pos[,2], pos[,2], ctB, ctA, w)
dy <- lmfity$coefficients[1]
dy
pos1[,1] <- pos1[,1]-dx
pos1[,2] <- pos1[,2]-dy

x <- c(pos1[,1], pos2[,1])
y <- c(pos1[,2], pos2[,2])
z <- c(rep(1, nrow(pos1)), rep(2, nrow(pos2)))
names(x) <- names(y) <- names(z) <- c(rownames(pos1), rownames(pos2))
plot3d(x,y,z, col=rainbow(2)[group], size=10)

pos <- rbind(pos1, pos2)
pos <- as.matrix(pos)
par(mfrow=c(1,1))
plotNetwork(pos, w, col=rainbow(2)[group], line.col = 'grey')
