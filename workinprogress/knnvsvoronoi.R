## Supp fig 2 KNN vs. voronoi

set.seed(0)
N <- 10^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
plotEmbedding(pos)

pos[,1] <- expPos(pos[,1])
pos[,2] <- expPos(pos[,2])

set.seed(0)
posj <- jitter(pos, amount=0.015)
plot(posj)

par(mfrow=c(1,3))

w <- getKnn(posj, 3)
plotNetwork(posj, w, main='KNN with K = 3')

w <- getKnn(posj, 6)
plotNetwork(posj, w, main='KNN with K = 6')

w <- voronoiAdjacency(posj)
plotNetwork(posj, w, main='Voronoi Adjacency')
