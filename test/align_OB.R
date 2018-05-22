# Sample analysis using OB data
library(MERingue)

# Data
data(cd)
data(pos)

cd <- as.matrix(cd[rownames(pos),])
sample_info <- pos
pos <- pos[,1:2]
hist(cd)

# split into left and right
plot(pos)
center = mean(pos[,1])
abline(v=center)
left <- pos[,1] < center
right <- pos[,1] > center

par(mfrow=c(1,2))
plot(pos[left,])
plot(pos[right,])

pos1 <- pos[left,]
pos2 <- pos[right,]
pos2[,1] <- -pos2[,1] # flip
cd1 <- cd[left,]
cd2 <- cd[right,]

# look for highly spatial genes
w <- getSpatialWeights(pos1, klist=3)
I <- getSpatialPatterns(t(cd1), w)
results <- I
vi <- results$p.adj < 0.05
vi[is.na(vi)] <- FALSE
table(vi)
results.sig1 <- rownames(results)[vi]
w <- getSpatialWeights(pos2, klist=3)
I <- getSpatialPatterns(t(cd2), w)
results <- I
vi <- results$p.adj < 0.05
vi[is.na(vi)] <- FALSE
table(vi)
results.sig2 <- rownames(results)[vi]

results.sig <- intersect(results.sig1, results.sig2)
gexp1 <- cd1[,results.sig]
gexp2 <- cd2[,results.sig]
par(mfrow=c(1,2))
plot(pos1, col=map2col(gexp1[,1]), pch=16)
plot(pos2, col=map2col(gexp2[,1]), pch=16)


# Plot
x <- c(pos1[,1], pos2[,1])
y <- c(pos1[,2], pos2[,2])
z <- c(rep(1, nrow(pos1)), rep(2, nrow(pos2)))
names(x) <- names(y) <- names(z) <- c(rownames(pos1), rownames(pos2))
#gexp <- c(gexp1, gexp2)
gexp <- rbind(gexp1, gexp2)
library(rgl)
plot3d(x,y,z, col=map2col(gexp[,1], colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)


# nearest neighbor in expression space
ctA <- rownames(pos1)
ctB <- rownames(pos2)
pos <- rbind(pos1, pos2)
pos <- as.matrix(pos)
group <- as.factor((c(ctA, ctB) %in% ctA)*1+1)
names(group) <- c(ctA, ctB)
table(group)
#gexpdf <- data.frame(c(gexp1, gexp2))
gexpdf <- data.frame(gexp)
#rownames(gexpdf) <- c(ctA, ctB)
w <- getMnn(ctA, ctB, pos=gexpdf, k=10)
plotNetwork(pos, w, col=rainbow(2)[group], line.col = 'grey')

plot3d(x,y,z, col=map2col(gexp[,1], colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)
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
plot3d(x,y,z, col=map2col(gexp[,1], colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)



# 2D interpolation
par(mfrow=c(2,2))
interpolate(pos1, gexp1[,1])
interpolate(pos2, gexp2[,1])
# 3D interpolation
library(deldir)
del <- deldir(x, y, z = z)
triangs <- do.call(rbind, triang.list(del))
plot3d(x,y,z, col=map2col(gexp[,1], colorRampPalette(c('blue', 'grey', 'red'))(100)), size=10)
triangles3d(triangs[, c("x", "y", "z")], col = "gray")
