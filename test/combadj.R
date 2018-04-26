# Test combining weights from multiple adjacency matrices

pos <- rbind(
  matrix(runif(1000, 0, 100), ncol=2),
  matrix(runif(1000, 100, 200), ncol=2)
)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- paste0('cell', 1:nrow(pos))

head(pos)

adj <- getAdj(pos, k=3)
plotNetwork(pos, adj)

# merge many adjs
adjList <- lapply(2:4, function(k) {
  adj <- getAdj(pos, k=k)
})
avgAdj <- Reduce("+", adjList) / length(adjList)

plotNetwork(pos, avgAdj, line.weight = 2)




cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

mat <- normalizeCounts(t(as.matrix(cd)),
                       verbose=FALSE)

results <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
  value <- mat[i,]
  ape::Moran.I(value, adj)
}, mc.cores=parallel::detectCores()-1))
rownames(results) <- rownames(mat)
results <- as.data.frame(results)
results$p.adj <- stats::p.adjust(results$p.value)

results2 <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
  value <- mat[i,]
  ape::Moran.I(value, avgAdj)
}, mc.cores=parallel::detectCores()-1))
rownames(results2) <- rownames(mat)
results2 <- as.data.frame(results2)
results2$p.adj <- stats::p.adjust(results2$p.value)
