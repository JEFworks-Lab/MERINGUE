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
pos <- jitter(pos_sub, factor = 5)

weight <- voronoiAdjacency(pos, plot=TRUE, filterDist=1)
rownames(weight) <- colnames(weight) <- rownames(pos)

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


# ## Simulate positive autocorrelation
# set.seed(3)
# gexp <- sort(rnorm(N*N))
# names(gexp) <- rownames(pos)[order(pos[,2])]
# gexp <- gexp[rownames(pos)]
# pv <- moranTest(gexp, weight)$p.value
# par(mfrow=c(1,1))
# plotNetwork(pos, weight, line.col='grey', main=pv, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
# #deldir::deldir(data$x,
# #               data$y,
# #               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey', main=pv)
# #box()
# points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)


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
points(pos, col=MERingue:::map2col(gexp), cex=2, pch=16)

lisaTest <- function (x, weight, na.rm=FALSE) {
  # unitization
  rs <- rowSums(weight)
  rs[rs == 0] <- 1
  weight <- weight/rs

  # calculate Ii
  n <- length(x)
  xx <- mean(x, na.rm = na.rm)
  z <- x - xx
  s2 <- sum(z^2, na.rm = na.rm)/n
  lz <- apply(weight, 1, function(w) sum(w*z))
  Ii <- (z/s2) * lz

  Wi <- rowSums(weight)
  E.Ii <- -Wi/(n - 1)

  b2 <- (sum(z^4, na.rm = na.rm)/n)/(s2^2)
  Wi2 <- apply(weight, 1, function(w) sum(w^2))
  A <- (n - b2)/(n - 1)
  B <- (2 * b2 - n)/((n - 1) * (n - 2))
  Sd.Ii <- sqrt(A * Wi2 + B * (Wi^2 - Wi2) - E.Ii^2)

  pv <- pnorm(Ii, mean = E.Ii, sd = Sd.Ii)
  pv <- 1 - pv

  return(data.frame(observed = Ii, expected = E.Ii, sd = Sd.Ii, p.value = pv))
}
lisa <- lisaTest(gexp, weight)
lisa
lp <- lisa$observed
#lp <- -log10(lisa$p.value)
#lp[lp > -log10(0.05)] <- -log10(0.05)
range(lp)
names(lp) <- rownames(lisa)
plotNetwork(pos, weight, line.col='grey', main='LISA', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(lp), cex=2, pch=16)

## Simulate negative autocorrelation
#set.seed(3)
pvs <- do.call(rbind, lapply(1:10000, function(i) {
  set.seed(i)
  gexp <- rnorm(N*N)
  names(gexp) <- rownames(pos)
  pv <- moranTest(gexp, weight)
  unlist(pv)
}))
rownames(pvs) <- 1:10000
pvs <- pvs[pvs[,1] < 0,]
pvs <- pvs[order(pvs[,1], decreasing=FALSE),]
head(pvs, n=10)

i <- rownames(pvs)[4]
set.seed(i)
gexp <- rnorm(N*N)
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
gexp <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp) <- rownames(pos)[order(pos[,2])]
gexp <- gexp[rownames(pos)]
par(mfrow=c(1,1))
plotNetwork(pos, weight, line.col='grey', main='Gene 1', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(pos, col=MERingue:::map2col(gexp), cex=1.5, pch=16)


set.seed(15)
cells1 <- paste0('cell', sample(1:N, N/2))
cells2 <- setdiff(paste0('cell', 1:N), cells1)

deldir::deldir(data$x,
               data$y,
               suppressMsge=TRUE, plotit = TRUE, sort=FALSE, col='grey')
points(pos[cells1,], col='red', cex=1.5, pch=16)
points(pos[cells2,], col='blue', cex=1.5, pch=16)

set.seed(1)
gexp1 <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
names(gexp1) <- c(cells1[order(pos[cells1,2])],
                  cells2[order(pos[cells2,2])])
#gexp1[cells2] <- 0
#gexp1 <- jitter(gexp1, factor=10)

gexp2 <- sort(log(abs(rnorm(N))+1), decreasing=FALSE)
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

moranTest(gexp1, weight)
moranTest(gexp2, weight)

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



par(mfrow=c(2,2))
hist(gexp1[cells1])
hist(gexp1[cells2])
hist(gexp2[cells1])
hist(gexp2[cells2])



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

