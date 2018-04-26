######################################### simulate spatially distributed data
pos <- rbind(
  matrix(runif(1000, 0, 100), ncol=2),
  matrix(runif(1000, 100, 200), ncol=2)
)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- paste0('cell', 1:nrow(pos))

head(pos)

########## gene expression is uniformly distributed (true negative)
set.seed(0)
value <- runif(nrow(pos))
names(value) <- rownames(pos)
plot(pos, col=map2col(value), pch=16)

# plot
adj <- getAdj(pos, k=3)
plotNetwork(pos, adj)

# calculate Moran's I
library(ape)
ape::Moran.I(value, adj)


########## gene expression is confounded with spatial pattern
set.seed(0)
value <- runif(nrow(pos))
value[pos$x > 100] <- value[pos$x > 100]*10
plot(pos, col=map2col(value), pch=16)

# plot
adj <- getAdj(pos, k=3)
plotNetwork(pos, adj)

# calculate Moran's I
library(ape)
ape::Moran.I(value, adj)


########## gene expression
set.seed(0)
value <- runif(nrow(pos))
value[pos$x > 150] <- value[pos$x > 150]*10
value[pos$x < 50] <- value[pos$x < 50]*10
plot(pos, col=map2col(value), pch=16)

# plot
adj <- getAdj(pos, k=3)
plotNetwork(pos, adj)

# calculate Moran's I
library(ape)
ape::Moran.I(value, adj)

