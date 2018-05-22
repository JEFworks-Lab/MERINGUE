# Simulation aligning different Z-planes

# Z1
set.seed(0)
N <- 100
pos1 <- cbind(rnorm(N), rnorm(N))
rownames(pos1) <- paste0('cell-z1-', 1:N)
colnames(pos1) <- c('x', 'y')
gexp1 <- pos1[,2]
plotEmbedding(pos1, colors=gexp1)

# Z2
set.seed(1)
N <- 200
pos2 <- cbind(rnorm(N), rnorm(N))
rownames(pos2) <- paste0('cell-z2-', 1:N)
colnames(pos2) <- c('x', 'y')
gexp2 <- pos2[,2]
plotEmbedding(pos2, colors=gexp2)

# Plot
library(scatterplot3d)
x <- c(pos1[,1], pos2[,1])
y <- c(pos1[,2], pos2[,2])
z <- c(rep(1, nrow(pos1)), rep(2, nrow(pos2)))
names(z) <- c(rownames(pos1), rownames(pos2))
gexp <- c(gexp1, gexp2)
scatterplot3d(x,y,z, pch=16, color=map2col(gexp, colorRampPalette(c('blue', 'grey', 'red'))(100)))
library(rgl)
plot3d(x,y,z, col=map2col(gexp, colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)

# Nearest neighbor across z-planes
ctA <- rownames(pos1)
ctB <- rownames(pos2)
pos <- rbind(pos1, pos2)
group <- as.factor((c(ctA, ctB) %in% ctA)*1+1)
names(group) <- c(ctA, ctB)
table(group)
w <- getMnn(ctA, ctB, pos, k=6)
plotNetwork(pos, w, col=rainbow(2)[group], line.col = 'grey')

# Visualize in 3D
library(rgl)
#plot3d(x=c(0,1),y=c(0,1),z=c(0,1), size=10)
#segments3d(
#  c(0,1),
#  c(0,1),
#  c(0,1)
#)
plot3d(x,y,z, col=map2col(gexp, colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)
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


# Look at expression correlation
lmfit <- plotNeighborBoxplot(gexp1, gexp2, ctA, ctB, w)
summary(lmfit)

# Look at spatial cross correlation
I = spatialCrossCor(c(gexp1, gexp2), c(gexp1, gexp2), ctA, ctB, w)
perm <- sapply(1:100, function(i) spatialCrossCor(sample(c(gexp1, gexp2)), sample(c(gexp1, gexp2)), ctA, ctB, w))
hist(perm, breaks=50, xlim=c(-5,5))
abline(v=I, col='red')
2 * pnorm(abs(I), mean(perm), sd(perm), lower.tail=FALSE) # two-sided p-value


#################### What happens if positions are shifted, how to optimally align?

# shift
pos2 <- pos2+10

x <- c(pos1[,1], pos2[,1])
y <- c(pos1[,2], pos2[,2])
z <- c(rep(1, nrow(pos1)), rep(2, nrow(pos2)))
names(z) <- c(rownames(pos1), rownames(pos2))
gexp <- c(gexp1, gexp2)
library(rgl)
plot3d(x,y,z, col=map2col(gexp, colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)

# nearest neighbor in expression space
ctA <- rownames(pos1)
ctB <- rownames(pos2)
pos <- rbind(pos1, pos2)
group <- as.factor((c(ctA, ctB) %in% ctA)*1+1)
names(group) <- c(ctA, ctB)
table(group)
gexpdf <- data.frame(c(gexp1, gexp2))
rownames(gexpdf) <- c(ctA, ctB)
w <- getMnn(ctA, ctB, pos=gexpdf, k=6)
plotNetwork(pos, w, col=rainbow(2)[group], line.col = 'grey')

plot3d(x,y,z, col=map2col(gexp, colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)
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

# look at positional differences
lmfitx <- plotNeighborBoxplot(pos[,1], pos[,1], ctA, ctB, w)
dx <- lmfitx$coefficients[1]

lmfity <- plotNeighborBoxplot(pos[,2], pos[,2], ctA, ctB, w)
dy <- lmfity$coefficients[1]

# align
pos2[,1] <- pos2[,1]-dx
pos2[,2] <- pos2[,2]-dy
x <- c(pos1[,1], pos2[,1])
y <- c(pos1[,2], pos2[,2])
z <- c(rep(1, nrow(pos1)), rep(2, nrow(pos2)))
names(z) <- c(rownames(pos1), rownames(pos2))
plot3d(x,y,z, col=map2col(gexp, colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)

