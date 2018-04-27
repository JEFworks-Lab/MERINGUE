
# Reanalysis 3

# mouse olfactory bulb
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# breast cancer
cd <- read.table('/Users/jefworks/Desktop/SpatialDE/Analysis/BreastCancer/data/Layer2_BC_count_matrix-1.tsv')
head(cd)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/BreastCancer/BC_sample_info.csv', row.names=1)
head(pos)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# seqFish
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/SeqFISH/exp_mat_43.csv', row.names=1, quote="\'", header=TRUE)
colnames(cd) <- gsub('X', '', colnames(cd))
cd <- t(cd)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/SeqFISH/sample_info_43.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]
head(cd)
head(pos)

# Clean
counts <- cleanCounts(t(cd), min.reads=10, min.detected=10)
mat <- normalizeCounts(counts)

# Analyze
adj <- getSpatialWeights(pos, plot=TRUE)
results <- getSpatialPatterns(mat, adj)
groups <- groupSigSpatialPatterns(pos, mat, results, k=5)

# cluster on interpolated?
vi <- results$p.adj < 0.05
results.sig <- rownames(results)[vi]
int <- do.call(cbind, lapply(results.sig, function(g) { as.vector(interpolate(pos, mat[g,], plot=FALSE)$z) }))
cv <- cor(as.matrix(int), use='pairwise.complete')
hc <- hclust(as.dist(1-cv))
groups <- as.factor(cutree(hc, 3))
par(mfrow=c(length(levels(groups)), 2), mar=rep(1,4))
prs <- lapply(levels(groups), function(g) {
  # summarize as first pc if more than 1 gene in group
  if(sum(groups==g)>1) {
    pc <- prcomp(mat[results.sig[groups==g],])
    pr <- pc$rotation[,1]
  } else {
    pr <- mat[results.sig[groups==g],]
  }
  interpolate(pos, pr, main=paste0("Pattern ", g, " : ", sum(groups==g), " genes"), plot=TRUE)
  return(pr)
})
names(prs) <- levels(groups)
# interpolation doesn't seem to help

######################################### simulate spatially distributed data
pos <- rbind(
  matrix(runif(1000, 0, 100), ncol=2),
  matrix(runif(1000, 100, 200), ncol=2)
)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- paste0('cell', 1:nrow(pos))

head(pos)

set.seed(0)
value1 <- runif(nrow(pos))
value2 <- runif(nrow(pos))
value2[pos$x > 100] <- value2[pos$x > 100]*10
value3 <- runif(nrow(pos))
value3[pos$x > 150] <- value3[pos$x > 150]*10
value3[pos$x < 50] <- value3[pos$x < 50]*10
value4 <- runif(nrow(pos))
value4[pos$x > 150] <- value4[pos$x > 150]*5
value4[pos$x < 50] <- value4[pos$x < 50]*5
mat <- rbind(value1, value2, value3, value4)
dim(mat)

# Analyze
adj <- getSpatialWeights(pos, klist=3, plot=TRUE)
results <- getSpatialPatterns(mat, adj)
results
groups <- groupSigSpatialPatterns(pos, mat, results, k=2, binSize=50)
