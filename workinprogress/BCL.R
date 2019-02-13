# Breast cancer layers analysis

library(MERingue)

data(BCL1)
data(BCL2)
data(BCL3)
data(BCL4)


## combine data
genes.have <- Reduce(intersect, list(
  rownames(BCL1$counts),
  rownames(BCL2$counts),
  rownames(BCL3$counts),
  rownames(BCL4$counts)

))
counts <- cbind(
  BCL1$counts[genes.have,],
  BCL2$counts[genes.have,],
  BCL3$counts[genes.have,],
  BCL4$counts[genes.have,]
)
colnames(counts) <- c(
  paste0('L1-', colnames(BCL1$counts)),
  paste0('L2-', colnames(BCL2$counts)),
  paste0('L3-', colnames(BCL3$counts)),
  paste0('L4-', colnames(BCL4$counts))
)
layer <-  c(
  rep('L1', ncol(BCL1$counts)),
  rep('L2', ncol(BCL2$counts)),
  rep('L3', ncol(BCL3$counts)),
  rep('L4', ncol(BCL4$counts))
)
names(layer) <- colnames(counts)
pos <- rbind(
  BCL1$pos[colnames(BCL1$counts),],
  BCL2$pos[colnames(BCL2$counts),],
  BCL3$pos[colnames(BCL3$counts),],
  BCL4$pos[colnames(BCL4$counts),]
)
rownames(pos) <- colnames(counts)

cc <- cleanCounts(counts, min.reads = 100, min.lib.size = 100)
mat <- normalizeCounts(cc, log=FALSE)
pos <- pos[colnames(mat),]
head(pos)

library(Matrix)
range(colSums(mat))
hist(colSums(as.matrix(mat)), breaks=20)
hist(colSums(as.matrix(cc)>0), breaks=20)

## batch correct between layers
#library(sva)
#mat.bc <- ComBat(as.matrix(mat), batch=factor(layer[colnames(mat)]))
#pos <- pos[colnames(mat.bc),]
mat.bc <- mat



par(mfrow=c(2,2), mar=rep(2,4))

library(jpeg)
pos1 <- BCL1$pos
rownames(pos1) <- paste0('L1-', rownames(pos1))
ima <- readJPEG("/Users/jefworks/Dropbox (HMS)/0Grad_School/Github/MERingue/data/HE_layer1_BC.jpg")
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(colnames(mat.bc), rownames(pos1)),], cex=1, lwd=3, col='black', pch=16)

library(jpeg)
pos1 <- BCL2$pos
rownames(pos1) <- paste0('L2-', rownames(pos1))
ima <- readJPEG("/Users/jefworks/Dropbox (HMS)/0Grad_School/Github/MERingue/data/HE_layer2_BC.jpg")
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(colnames(mat.bc), rownames(pos1)),], cex=1, lwd=3, col='black', pch=16)

library(jpeg)
pos1 <- BCL3$pos
rownames(pos1) <- paste0('L3-', rownames(pos1))
ima <- readJPEG("/Users/jefworks/Dropbox (HMS)/0Grad_School/Github/MERingue/data/HE_layer3_BC.jpg")
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(colnames(mat.bc), rownames(pos1)),], cex=1, lwd=3, col='black', pch=16)

library(jpeg)
pos1 <- BCL4$pos
rownames(pos1) <- paste0('L4-', rownames(pos1))
ima <- readJPEG("/Users/jefworks/Dropbox (HMS)/0Grad_School/Github/MERingue/data/HE_layer4_BC.jpg")
plotEmbedding(pos1)
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
points(pos1[intersect(colnames(mat.bc), rownames(pos1)),], cex=1, lwd=3, col='black', pch=16)



##################### Get highly spatially variable genes for each layer individually
par(mfrow=c(1,1))
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
w <- voronoiAdjacency(pos[subset,], plot=TRUE)
plotNetwork(pos=pos[subset,], adj=w)
I <- getSpatialPatterns(mat.bc[, subset], w)
results.filter <- filterSpatialPatterns(mat = mat.bc[, subset],
                                        I = I,
                                        w=w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
scc <- spatialCrossCorMatrix(mat.bc[results.filter, subset], w)
ggroup <- groupSigSpatialPatterns(pos=pos[subset,],
                                  mat=as.matrix(mat.bc[results.filter, subset]),
                                  scc=scc,
                                  plot=TRUE,
                                  hclustMethod='ward.D',
                                  power=1,
                                  deepSplit=0)

## plot one gene
g <- names(ggroup$groups)[ggroup$groups==levels(ggroup$groups)[1]][2]
g
invisible(interpolate(pos[subset,],
                      mat.bc[g, subset],
                      main = g))
invisible(interpolate(BCL2$pos,
                      BCL2$counts[g, ],
                      main = g))

####### repeat for all layers
helper <- function(subset) {
  w <- voronoiAdjacency(pos[subset,])
  plotNetwork(pos=pos[subset,], adj=w)
  I <- getSpatialPatterns(mat.bc[, subset], w)

  out <- filterSpatialPatterns(mat = mat.bc[, subset],
                                          I = I,
                                          w = w,
                                          adjustPv = TRUE,
                                          alpha = 0.05,
                                          minPercentCells = 0.05,
                                          verbose = TRUE,
                                          details=TRUE)

  results.filter <- filterSpatialPatterns(mat = mat.bc[, subset],
                                          I = I,
                                          w = w,
                                          adjustPv = TRUE,
                                          alpha = 0.05,
                                          minPercentCells = 0.05,
                                          verbose = TRUE)
  scc <- spatialCrossCorMatrix(mat.bc[results.filter, subset], w)
  ggroup <- groupSigSpatialPatterns(pos=pos[subset,],
                                    mat=as.matrix(mat.bc[results.filter, subset]),
                                    scc=scc,
                                    plot=TRUE,
                                    hclustMethod='ward.D',
                                    power=1,
                                    deepSplit=0)
  list(I=I, sig.genes=results.filter, ggroup=ggroup, out=out)
}

subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
L1 <- helper(subset)

subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
L2 <- helper(subset)

subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
L3 <- helper(subset)

subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
L4 <- helper(subset)

## compare
library(UpSetR)
sig <- list('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes)
upset(fromList(sig))

length(unique(unlist(sig)))

pv <- cbind(L1$I$p.adj, L2$I$p.adj, L3$I$p.adj, L4$I$p.adj)
colnames(pv) <- c('L1', 'L2', 'L3', 'L4')
rownames(pv) <- rownames(mat.bc)
head(pv)

pv[Reduce(intersect, sig),]

head(L1$out)
write.csv(L1$out, file="L1.csv")
write.csv(L2$out, file="L2.csv")
write.csv(L3$out, file="L3.csv")
write.csv(L4$out, file="L4.csv")


################### Use 3D information
pos3d <- cbind(pos, layer = as.numeric(factor(layer[rownames(pos)])))
head(pos3d)
pos3d <- pos3d[colnames(mat.bc),]
dim(pos3d)

plot(BCL1$pos)
plot(BCL2$pos)
plot(BCL3$pos)
plot(BCL4$pos)

library(scatterplot3d)
par(mfrow=c(1,1))
scatterplot3d(pos3d, color=rainbow(4)[as.factor(pos3d[,3])])

library(rgl)
plot3d(pos3d[,1],pos3d[,2],pos3d[,3], col=rainbow(4)[as.factor(pos3d[,3])], size=10, xlab='', ylab='', zlab='')
planes3d(0, 0, 1, -1, alpha = 0.1)
planes3d(0, 0, 1, -2, alpha = 0.1)
planes3d(0, 0, 1, -3, alpha = 0.1)
planes3d(0, 0, 1, -4, alpha = 0.1)


between <- lapply(sort(unique(pos3d[,3]))[-1], function(i) {
  ctA <- rownames(pos3d)[pos3d[,3]==(i-1)]
  ctB <- rownames(pos3d)[pos3d[,3]==i]
  if(length(ctA)==0 | length(ctB)==0) {
    wa <- NA
  } else {
    wa <- getMnn(ctA, ctB, pos=pos3d[c(ctA,ctB),1:2], k=3)
  }
  return(wa)
})
w <- matrix(0, nrow=nrow(pos3d), ncol=nrow(pos3d))
rownames(w) <- colnames(w) <- rownames(pos3d)
## across laters
invisible(lapply(between, function(wa) {
  w[rownames(wa), colnames(wa)] <<- wa
}))


# 3d plot test
x <- pos3d[,1]
y <- pos3d[,2]
z <- pos3d[,3]
names(x) <- names(y) <- names(z) <- rownames(pos3d)
plot3d(-x, y, z, col=rainbow(4)[as.factor(pos3d[,3])], size=15, alpha = 0.50,
       xlab='', ylab='', zlab='', axes=FALSE, expand = 1.1)
rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('white'), expand=1.1)

idx <- which(w>0, arr.ind = T)
for(i in sample(seq_len(nrow(idx)), 500)) {
  n1 <- rownames(w)[idx[i,1]]
  n2 <- rownames(w)[idx[i,2]]
  segments3d(
    x=c(-x[n1], -x[n2]),
    y=c(y[n1], y[n2]),
    z=c(z[n1], z[n2]),
    col='lightgrey'
  )
}

# 3d plot test 2
plt <- scatterplot3d(-pos3d[,3],pos3d[,1],pos3d[,2], color=rainbow(4)[as.factor(pos3d[,3])], pch=16, zlim=c(-15,15), ylim=c(-15,15), axis=FALSE)
vi <- which(w==1, arr.ind=TRUE)
from <- rownames(w)[vi[,'row']]
to <- colnames(w)[vi[,'col']]
head(from)
head(to)
sapply(sample(seq_along(from), 500), function(i) {
  plt$points3d(x=c(-pos3d[from[i],3], -pos3d[to[i],3]),
               y=c(pos3d[from[i],1], pos3d[to[i],1]),
               z=c(pos3d[from[i],2], pos3d[to[i],2]),
               type="l", col=rgb(0.8,0.8,0.8,0.4), lwd=2)
})


getCrossLayerNeighbors <- function(layers, k=3) {
  subset <- unlist(lapply(1:length(layers), function(i) {
    paste0(i, rownames(layers[[i]]))
  }))
  between <- lapply(1:(length(layers)-1), function(i) {
    pos1 <- layers[[i]]
    pos2 <- layers[[i+1]]
    ctA <- paste0(i, rownames(pos1))
    rownames(pos1) <- ctA
    ctB <- paste0(i+1, rownames(pos2))
    rownames(pos2) <- ctB
    pos <- rbind(pos1, pos2)
    if(length(ctA)==0 | length(ctB)==0) {
      wa1 <- NA
    } else {
      pos <- as.matrix(pos)
      wa1 <- getMnn(ctA, ctB, pos=pos, k=k)
    }
    return(wa1)
  })

  w <- matrix(0, nrow=length(subset), ncol=length(subset))
  rownames(w) <- colnames(w) <- subset
  ## across layers
  invisible(lapply(between, function(wa1) {
    w[rownames(wa1), colnames(wa1)] <<- wa1
  }))

  return(w)
}

layers <- list(pos[colnames(mat.bc)[grepl('L1', colnames(mat.bc))],],
               pos[colnames(mat.bc)[grepl('L2', colnames(mat.bc))],],
               pos[colnames(mat.bc)[grepl('L3', colnames(mat.bc))],],
               pos[colnames(mat.bc)[grepl('L4', colnames(mat.bc))],]
)
w <- getCrossLayerNeighbors(layers, k=6)
I <- getSpatialPatterns(mat.bc, w)
out <- filterSpatialPatterns(mat = mat.bc,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0, # just confirming
                                        verbose = TRUE, details=TRUE)
write.csv(out, file="cross.csv")
results.filter <- rownames(out)
length(results.filter)
cross <- list(I=I, sig.genes=results.filter, out=out)

## compare
library(UpSetR)
sig <- list('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes, 'Cross'=cross$sig.genes)
?upset
upset(fromList(sig), keep.order = TRUE, order.by = 'degree')
lapply(sig, length)

pv <- cbind(L1$I$p.adj, L2$I$p.adj, L3$I$p.adj, L4$I$p.adj, cross$I$p.adj)
colnames(pv) <- c('L1', 'L2', 'L3', 'L4', 'Cross')
rownames(pv) <- rownames(mat.bc)
head(pv)

shared <- pv[Reduce(intersect, sig),]
shared


subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
w <- voronoiAdjacency(pos[subset,])
scc <- spatialCrossCorMatrix(mat.bc[rownames(shared), subset], w)
ggroup <- groupSigSpatialPatterns(pos=pos[subset,],
                                  mat=as.matrix(mat.bc[rownames(shared), subset]),
                                  scc=scc,
                                  plot=TRUE,
                                  hclustMethod='ward.D',
                                  power=1,
                                  deepSplit=0)

par(mfrow=c(2,2))
#g <- rownames(shared)[16]
intersect(rownames(shared), names(ggroup$groups)[which(ggroup$groups==3)])
#intersect(consistent.genes, names(ggroup$groups)[which(ggroup$groups==2)])
#g <-  "PRSS23"
g <- "RPL13"
zlim <- c(-0.5,0.5)
subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)


############## Look for genes that are consistent across at least 2 layers
consistent.genes <- intersect(Reduce(union, list('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes)), cross$sig.genes)
length(consistent.genes)

#bar <- setdiff(intersect(L1$sig.genes, L3$sig.genes), consistent.genes)
bar <- setdiff(intersect(L3$sig.genes, L4$sig.genes), consistent.genes)
bar
par(mfrow=c(2,2))
g <- bar[1]
-log10(foo[g,])
zlim <- c(-1,1)
subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)




## genes from just looking at individual slices
foo <- unlist(c('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes))
foo <- foo[duplicated(foo)]
length(foo)
setdiff(foo, consistent.genes) ## difference: genes that are clustered in different slices but not in a coordinated manner

inconsistent.genes <- setdiff(Reduce(union, list('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes)), cross$sig.genes)
foo <- pv[inconsistent.genes,]
vi <- apply(foo, 1, function(x) sum(x<0.05)>1) # remove those unique to a layer
table(vi)
bar <- foo[vi,]
-log10(bar[1,])

par(mfrow=c(2,2))
#g <- rownames(bar)[3]
g <- "COL12A1"
zlim <- c(-0.5,0.5)
subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)


############## Look at a specific gene
g <-  "RPL29"
par(mfrow=c(1,1))
#scatterplot3d(pos3d, color=map2col(winsorize(log10(mat[g,]+1), 0.05)), pch=16)
#scatterplot3d(pos3d[,3],pos3d[,2],pos3d[,1], color=MERingue:::map2col(winsorize(mat.bc[g,], 0.05)), pch=16, zlim=c(-15,15), ylim=c(-15,15), axis=FALSE)
scatterplot3d(-pos3d[,3],pos3d[,1],pos3d[,2], color=MERingue:::map2col(winsorize(mat.bc[g,], 0.05)), pch=16, zlim=c(-15,15), ylim=c(-15,15), axis=FALSE)

pv[g,]

## genes from just looking at individual slices
foo <- unique(unlist(c('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes)))
length(foo)
setdiff(cross$sig.genes, foo) # things we gain from cross layer analysis not apparent in intersection
g <-"ACACA"
par(mfrow=c(1,1))
scatterplot3d(pos3d[,1],pos3d[,2],pos3d[,3], color=MERingue:::map2col(winsorize(mat.bc[g,])), pch=16, xlim=c(-15,15), ylim=c(-15,15), axis=FALSE)

foo <- unlist(c('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes))
foo <- foo[duplicated(foo)]
setdiff(foo, cross$sig.genes)
g <-"ACACA"
L2$I[g,]
cross$I[g,]
par(mfrow=c(1,1))
scatterplot3d(pos3d[,1],pos3d[,2],pos3d[,3], color=MERingue:::map2col(winsorize(mat.bc[g,])), pch=16, xlim=c(-15,15), ylim=c(-15,15), axis=FALSE)


par(mfrow=c(2,2), mar=rep(1,4))
g <- "COX6B1"
L1$I[g,]
L2$I[g,]
L3$I[g,]
L4$I[g,]
cross$I[g,]
subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g)
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g)
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g)
subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g)



##### consistent
#gs <- rownames(pv[Reduce(intersect, sig),])
gs <- consistent.genes
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
w <- voronoiAdjacency(pos[subset,], filterDist=1.75)
scc <- spatialCrossCorMatrix(mat.bc[gs, subset], w)
ggroup <- groupSigSpatialPatterns(pos=pos[subset,],
                                  mat=as.matrix(mat.bc[gs, subset]),
                                  scc=scc,
                                  plot=TRUE,
                                  hclustMethod='ward.D',
                                  power=1,
                                  deepSplit=0)



## unique to a slice
#gs <- setdiff(L1$sig.genes, Reduce(union, list('L4'=L4$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'cross'=cross$sig.genes)))
#gs <- setdiff(L2$sig.genes, Reduce(union, list('L4'=L4$sig.genes, 'L1'=L1$sig.genes, 'L3'=L3$sig.genes, 'cross'=cross$sig.genes)))
#gs <- setdiff(L4$sig.genes, Reduce(union, list('L2'=L2$sig.genes, 'L1'=L1$sig.genes, 'L3'=L3$sig.genes, 'cross'=cross$sig.genes)))
gs <- setdiff(L3$sig.genes, Reduce(union, list('L2'=L2$sig.genes, 'L1'=L1$sig.genes, 'L4'=L4$sig.genes, 'cross'=cross$sig.genes)))
gs
#gs <- names(L1$ggroup$groups)[which(L1$ggroup$groups == 1)]
#gs <- intersect(Reduce(intersect, list('L1'=L1$sig.genes, 'L2'=L2$sig.genes, 'L3'=L3$sig.genes, 'L4'=L4$sig.genes)), cross$sig.genes)
#gs
#gs <- setdiff(intersect(L1$sig.genes, L2$sig.genes), cross$sig.genes)
#gs

#gs <- c(consistent.genes, inconsistent.genes)
#subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
#subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
#subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
scc <- spatialCrossCorMatrix(mat.bc[gs, subset], voronoiAdjacency(pos[subset,]))
ggroup <- groupSigSpatialPatterns(pos=pos[subset,],
                                  mat=t(scale(t(as.matrix(mat.bc[gs, subset])))),
                                  scc=scc,
                                  plot=TRUE,
                                  hclustMethod='ward.D',
                                  power=1,
                                  deepSplit=0)


gs <- names(ggroup$groups)[ggroup$groups==4]

par(mfrow=c(2,2))
g <- gs[4]
foo[g,]
zlim <- c(-0.5,0.5)
subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)
subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=scale(mat[g,subset])[,1], cex=2, main=g, zlim=zlim, alpha=1)






gs1 <- intersect(gs, names(ggroup$groups)[which(ggroup$groups == 4)])
gs1

lapply(gs1, function(g) {


par(mfrow=c(2,2))
#g <- gs1[1]
#gexp <- scale(log10(mat[g,]+1))[,1]
gexp <- scale(mat[g,])[,1]
#gexp <- mat[g,]
range(gexp)
gexp[gexp > 0.4] <- 0.4
gexp[gexp < -0.4] <- -0.4

subset <- colnames(mat.bc)[grepl('L1', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=gexp[subset], main=g, alpha=1)
subset <- colnames(mat.bc)[grepl('L2', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=gexp[subset],  main=g, alpha=1)
subset <- colnames(mat.bc)[grepl('L3', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=gexp[subset], main=g, alpha=1)
subset <- colnames(mat.bc)[grepl('L4', colnames(mat.bc))]
plotEmbedding(pos[subset,], col=gexp[subset], main=g, alpha=1)

})
