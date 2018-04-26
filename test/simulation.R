# simulate spatially distributed data
pos <- rbind(
  matrix(rnorm(100, 0, 10), ncol=2),
  matrix(rnorm(100, 100, 10), ncol=2)
)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- paste0('cell', 1:nrow(pos))

head(pos)

# gene expression is uniformly distributed
set.seed(0)
value <- rnorm(nrow(pos))
names(value) <- rownames(pos)

# plot
adj1 <- getAdj(pos, k=3)
adj2 <- getAdj(pos, k=10)
plotNetwork(pos, adj1)
plotNetwork(pos, adj2)

# calculate Moran's I
library(ape)
ape::Moran.I(value, adj1)
ape::Moran.I(value, adj2)

plot(pos, col=map2col(value), pch=16)
