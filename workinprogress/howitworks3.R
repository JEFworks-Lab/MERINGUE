# Make Fig 1A

par(mfrow=c(1,1))

set.seed(0)
N <- 10^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
pos <- pos*1.5
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
plotEmbedding(pos)

# induce non-homogeneity
pos_sub <- pos
pos_sub[,1] <- 1.1^abs(pos_sub[,1])
pos_sub[,2] <- 1.1^abs(pos_sub[,2])
plotEmbedding(pos_sub)
pos <- pos_sub

weight <- voronoiAdjacency(pos, plot=TRUE)

makeixy <- function(data, formula, scale){
  m = model.frame(formula, data=data)
  if(ncol(m)!=3){
    stop("incorrect adjacency formula: id~x+y needed")
  }
  names(m)=c("id","x","y")
  m[,2]=m[,2]/scale
  m[,3]=m[,3]/scale
  m
}
data <- data.frame('i'=1:nrow(pos), pos)
data <- makeixy(data, formula=i~x+y, scale=1)
dd1 = deldir::deldir(data$x,
                     data$y,
                     suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey')
box()
adj = weight
idx <- which(adj>0, arr.ind = T)
#line.col= rgb(0.5,0.5,0.5,0.5)
line.col= 'red'
line.power = 2
for(i in seq_len(nrow(idx))) {
  lines(
    c(pos[idx[i,1],1], pos[idx[i,2],1]),
    c(pos[idx[i,1],2], pos[idx[i,2],2]),
    col=line.col,
    lwd=adj[idx]*line.power,
  )
}
points(pos, pch=16, cex=2)



## Simulate correlation
par(mfrow=c(1,1))
set.seed(1)
cells1 <- paste0('cell', sample(1:N, N/2))
cells2 <- setdiff(paste0('cell', 1:N), cells1)

deldir::deldir(data$x,
               data$y,
               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey')
points(pos[cells1,], col='darkorange', cex=1.5, pch=16)
points(pos[cells2,], col='darkgreen', cex=1.5, pch=16)

set.seed(1)
gexp <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp) <- rownames(pos)[order(pos[,2])]
gexp <- gexp[rownames(pos)]

gexp1 <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp1) <- c(cells1[order(pos[cells1,2])],
                  cells2[order(pos[cells2,2])])

gexp2 <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp2) <- c(cells2[order(pos[cells2,2])],
                  cells1[order(pos[cells1,2])])
gexp1 <- gexp1[rownames(pos)]
gexp2 <- gexp2[rownames(pos)]

par(mfrow=c(1,1), mar=rep(5,4))
#plotNetwork(pos, weight, line.col='grey', main='Gene 0', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#points(pos, col=MERingue:::map2col(gexp), cex=2, pch=16)
plotNetwork(pos, weight, line.col='grey', main='Gene 1', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(gexp1), cex=2, pch=16)
plotNetwork(pos, weight, line.col='grey', main='Gene 2', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(gexp2), cex=2, pch=16)

moranTest(gexp, weight)
moranTest(gexp1, weight)
moranTest(gexp2, weight)

mat1 <- rbind(gexp, gexp1, gexp2)

## Simulate correlation 2
set.seed(15)
gexp <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp) <- rownames(pos)[order(pos[,1])]
gexp <- gexp[rownames(pos)]

gexp1 <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp1) <- c(cells1[order(pos[cells1,1])],
                  cells2[order(pos[cells2,1])])

gexp2 <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp2) <- c(cells2[order(pos[cells2,1])],
                  cells1[order(pos[cells1,1])])

gexp1 <- gexp1[rownames(pos)]
gexp2 <- gexp2[rownames(pos)]

plotNetwork(pos, weight, line.col='grey', main='Gene 0', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(gexp), cex=2, pch=16)
plotNetwork(pos, weight, line.col='grey', main='Gene 1', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(gexp1), cex=2, pch=16)
plotNetwork(pos, weight, line.col='grey', main='Gene 2', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(gexp2), cex=2, pch=16)

moranTest(gexp, weight)
moranTest(gexp1, weight)
moranTest(gexp2, weight)

mat2 <- rbind(gexp, gexp1, gexp2)
rownames(mat2) <- paste0('pattern 2', rownames(mat2))

mat <- rbind(mat1, mat2)
scc <- spatialCrossCorMatrix(mat, weight)
heatmap(scc, col=colorRampPalette(c('white', 'black'))(100))
#scc <- dist(mat)
groupSigSpatialPatterns(pos, mat, scc)


gexp1 <- mat1[2,]
gexp2 <- mat1[3,]
summary(lm(gexp1~gexp2))
plot(gexp1, gexp2, pch=".")
points(gexp1[cells1], gexp2[cells1], pch=15, cex=1.5, col='darkgreen')
points(gexp1[cells2], gexp2[cells2], pch=17, cex=1.5, col='darkorange')
abline(lm(gexp1~gexp2), col='red')
text(0.5,0.5, paste0('cor ', cor.test(gexp1, gexp2)$estimate), col='red')

I <- spatialIntraCrossCor(gexp1, gexp2, weight=weight)
I

# random permutation
perm <- sapply(1:1000, function(i) spatialIntraCrossCor(sample(gexp1), sample(gexp2), weight))
hist(perm, breaks=50, xlim=c(-1,1))
abline(v=I, col='red')
pv <- 2 * pnorm(abs(I), mean(perm), sd(perm), lower.tail=FALSE) # two-sided p-value
text(paste0('p-value ', pv), col='red', x=0, y=50)

plotNeighborBoxplot(gexp1, gexp2, cells1, cells2, weight)
gexp1
gexp2

I <- spatialCrossCor(gexp1, gexp2, cells1, cells2, weight=weight)
I






## Simulate negative autocorrelation
#set.seed(3)
pvs <- do.call(rbind, lapply(1:100000, function(i) {
  set.seed(i)
  gexp <- rnorm(N)
  names(gexp) <- rownames(pos)
  pv <- moranTest(gexp, weight)
  unlist(pv)
}))
rownames(pvs) <- 1:100000
pvs <- pvs[pvs[,1] < 0,]
pvs <- pvs[order(pvs[,1], decreasing=FALSE),]
head(pvs, n=10)

par(mfrow=c(2,3))
lapply(1:6, function(j) {
  i <- rownames(pvs)[j]
  set.seed(i)
  gexp <- rnorm(N)
  #hist(gexp)
  names(gexp) <- rownames(pos)
  pv <- moranTest(gexp, weight)$observed
  plotNetwork(pos, weight, line.col='grey', main=pv, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
  points(pos, col=MERingue:::map2col(gexp), cex=2, pch=16)
  #plotEmbedding(pos, col=gexp, cex=2)
})




gexpA <- gexp1
gexpB <- gexp2
gexpB[cells1] <- 0
#gexpB[cells2] <- scale(gexpB[cells2])[,1]
gexpA[cells2] <- 0
#gexpA[cells1] <- scale(gexpA[cells1])[,1]

w <- getMnn(cells1, cells2, pos, k=5)

par(mfrow=c(3,3), mar=rep(5,4))
plotNetwork(pos, w, line.col='lightgrey', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos[cells1,], col='darkgreen', cex=2, pch=16)
points(pos[cells2,], col='darkorange', cex=2, pch=15)
spatialCrossCor(gexp1, gexp2, cells1, cells2, w)/2

plot(gexpA, gexpB, pch='.',  xlab='Gene A', ylab='Gene B', main='Correlation')
points(gexpA[cells1], gexpB[cells1], col='darkgreen', cex=2, pch=16)
points(gexpA[cells2], gexpB[cells2], col='darkorange', cex=2, pch=15)
abline(lm(gexpB~gexpA), col='red')
text(0.3,0.3, paste0('R = ', cor(gexpA, gexpB)), col='red')

#plot(pos, xlab=NA, ylab=NA, xaxt='n', yaxt='n', pch=".", main='Gene 1')
plot(pos[cells1,], col=MERingue:::map2col(gexpA)[cells1], cex=2, pch=16, xlab=NA, ylab=NA, xaxt='n', yaxt='n',
     main = 'Gene A in Cell-Type A')
plot(pos[cells2,], col=MERingue:::map2col(gexpB)[cells2], cex=2, pch=15, xlab=NA, ylab=NA, xaxt='n', yaxt='n',
     main = 'Gene B in Cell-Type B')

#plot(pos, xlab=NA, ylab=NA, xaxt='n', yaxt='n', pch=".", main='Gene 2')
plot(pos[cells1,], col=MERingue:::map2col(gexpB)[cells1], cex=2, pch=16, xlab=NA, ylab=NA, xaxt='n', yaxt='n',
     main = 'Gene B in Cell-Type A')
plot(pos[cells2,], col=MERingue:::map2col(gexpA)[cells2], cex=2, pch=15, xlab=NA, ylab=NA, xaxt='n', yaxt='n',
     main = 'Gene A in Cell-Type B')

plotNeighborBoxplot(gexp1, gexp2, cells1, cells2, w)

cor(gexpA, gexpB)
spatialCrossCor(gexp1, gexp2, cells1, cells2, w)/2

bg <- unlist(lapply(1:1000, function(i) {
  set.seed(i)
  bg <- spatialCrossCor(sample(gexp1, length(gexp1)), sample(gexp2, length(gexp2)), cells1, cells2, w)/2
}))
hist(bg, breaks=20, xlim=c(-1, 1), xlab='iSCC')
abline(v = spatialCrossCor(gexp1, gexp2, cells1, cells2, w)/2, col='red')


spatialIntraCrossCor(gexp1, gexp2, weight)/2
bg <- unlist(lapply(1:1000, function(i) {
  set.seed(i)
  bg <- spatialIntraCrossCor(sample(gexp1, length(gexp1)), sample(gexp2, length(gexp2)), weight)/2
}))
hist(bg, breaks=20, xlim=c(-1, 1), xlab='iSCC')
abline(v = spatialCrossCor(gexp1, gexp2, cells1, cells2, weight)/2, col='red')

rownames(weight) <- colnames(weight) <- c(cells1, cells2)
plotNeighborBoxplot(gexp1, gexp2, c(cells1, cells2), c(cells1, cells2), weight)
