# Make Fig 1A

set.seed(2)
N <- 100
pos <- cbind(abs(rnorm(N, 0, 1))^1.1, rnorm(N, 0, 1))
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
pos
weight <- voronoiAdjacency(pos, plot=TRUE, filterDist=1)

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
line.col= rgb(0.5,0.5,0.5,0.5)
line.power = 2
for(i in seq_len(nrow(idx))) {
  lines(
    c(pos[idx[i,1],1], pos[idx[i,2],1]),
    c(pos[idx[i,1],2], pos[idx[i,2],2]),
    col=line.col,
    lwd=adj[idx]*line.power,
  )
}
points(pos, pch=16)


## Simulate positive autocorrelation
set.seed(3)
gexp <- sort(rnorm(N))
names(gexp) <- rownames(pos)[order(pos[,2])]
gexp <- gexp[rownames(pos)]
pv <- moranTest(gexp, weight)$p.value
par(mfrow=c(1,1))
plotNetwork(pos, weight, line.col='grey', main=pv, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#deldir::deldir(data$x,
#               data$y,
#               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', main=pv)
#box()
points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)


## Simulate positive autocorrelation
set.seed(3)
gexp <- sort(rnorm(N))
names(gexp) <- rownames(pos)[order((pos[,2]-mean(pos[,2]))*(pos[,1]-mean(pos[,1])))]
gexp <- gexp[rownames(pos)]
pv <- moranTest(gexp, weight)$p.value
par(mfrow=c(1,1))
plotNetwork(pos, weight, line.col='grey', main=pv, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#deldir::deldir(data$x,
#               data$y,
#               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', main=pv)
#box()
points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)



## Simulate negative autocorrelation
#set.seed(3)
pvs <- do.call(rbind, lapply(1:10000, function(i) {
  set.seed(i)
  gexp <- rnorm(N)
  names(gexp) <- rownames(pos)
  pv <- moranTest(gexp, weight)
  unlist(pv)
}))
rownames(pvs) <- 1:10000
pvs <- pvs[pvs[,1] < 0,]
pvs <- pvs[order(pvs[,1], decreasing=FALSE),]
head(pvs, n=10)

i <- rownames(pvs)[1]
set.seed(i)
gexp <- rnorm(N)
hist(gexp)
names(gexp) <- rownames(pos)
pv <- moranTest(gexp, weight)$observed

par(mfrow=c(1,1))
plotNetwork(pos, weight, line.col='grey', main=pv, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#deldir::deldir(data$x,
#               data$y,
#               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', main=pv)
points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)





## Simulate cross correlation
set.seed(15)
cells1 <- paste0('cell', sample(1:N, N/2))
cells2 <- setdiff(paste0('cell', 1:N), cells1)

deldir::deldir(data$x,
               data$y,
               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey')
points(pos[cells1,], col='red', cex=1.5, pch=16)
points(pos[cells2,], col='blue', cex=1.5, pch=16)

set.seed(0)
gexp1 <- sort(log(abs(rnorm(N))+1), decreasing=TRUE)
names(gexp1) <- c(cells1[order(pos[cells1,2])],
                  cells2[order(pos[cells2,2])])
#gexp1[cells2] <- 0
#gexp1 <- jitter(gexp1, factor=10)

gexp2 <- sort(log(abs(rnorm(N))+1), decreasing=TRUE)
names(gexp2) <- c(cells2[order(pos[cells2,2])],
                  cells1[order(pos[cells1,2])])
#gexp2[cells1] <- 0
#gexp2 <- jitter(gexp2, factor=10)

gexp1 <- gexp1[rownames(pos)]
gexp2 <- gexp2[rownames(pos)]

par(mfrow=c(1,1))
plotNetwork(pos, weight, line.col='grey', main='Gene 1', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#deldir::deldir(data$x,
#               data$y,
#               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', main='Gene 1')
#points(pos[cells1,], col=MERingue:::map2col(gexp1[cells1]), cex=1.5, pch=16)
points(pos, col=MERingue:::map2col(gexp1), cex=1.5, pch=16)
#box()

plotNetwork(pos, weight, line.col='grey', main='Gene 2', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#deldir::deldir(data$x,
#               data$y,
#               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', main='Gene 2')
#points(pos[cells2,], col=MERingue:::map2col(gexp2[cells2]), cex=1.5, pch=16)
points(pos, col=MERingue:::map2col(gexp2), cex=1.5, pch=16)
#box()


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

I <- spatialCrossCor(gexp1, gexp2, cells1, cells2, weight=weight)
I




hist(gexp1[cells1])
hist(gexp1[cells2])



set.seed(11)
gexp <- sort(rnorm(N))
names(gexp) <- rownames(pos)[order(pos[,2]-pos[,1])]
gexp <- gexp[rownames(pos)]
#gexp[pos[,1] < 1] <- 0
pv <- moranTest(gexp, weight)$p.value
par(mfrow=c(1,1))
deldir::deldir(data$x,
               data$y,
               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', frame.plot=TRUE)
points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)
gexp1 <- scale(gexp)[,1]

set.seed(7)
gexp <- sort(rnorm(N))
names(gexp) <- rownames(pos)[order(pos[,2]+pos[,1])]
gexp <- gexp[rownames(pos)]
#gexp[pos[,1] < 1] <- 0
pv <- moranTest(gexp, weight)$p.value
par(mfrow=c(1,1))
deldir::deldir(data$x,
               data$y,
               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', frame.plot=TRUE)
points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)
gexp2 <- scale(gexp)[,1]

plot(gexp1, gexp2)
cor.test(gexp1, gexp2)

spatialIntraCrossCor(gexp1, gexp2, weight=weight)






plotNetwork(pos, weight, line.col='grey', main='Cell Types', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos[cells1,], col='darkgreen', cex=1.5, pch=15)
points(pos[cells2,], col='darkorange', cex=1.5, pch=17)

