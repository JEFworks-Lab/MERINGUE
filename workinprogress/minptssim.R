## Minimum number of points simulation
library(MERingue)
K <- 20

set.seed(0)
N <- K^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
pos <- pos*1.5
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
plotEmbedding(pos)

weight <- voronoiAdjacency(pos, plot=FALSE)

## random gene expression
M <- 10000
rand <- matrix(rnorm(N*M), ncol=N)
dim(rand)
colnames(rand) <- rownames(pos)
rownames(rand) <- paste0('gene', 1:M)

## test
getSpatialPatterns <- function(mat, adj, verbose=TRUE) {

  # Calculate Moran's I for each gene
  results <- MERingue:::getSpatialPatterns_C(as.matrix(mat), as.matrix(adj), verbose)
  colnames(results) <- c('observed', 'expected', 'sd')
  rownames(results) <- rownames(mat)
  results <- as.data.frame(results)

  # Get p-value
  # always assume we want greater autocorrelation
  pv <- 1 - pnorm(results$observed, mean = results$expected, sd = results$sd)
  results$p.value <- pv;
  # multiple testing correction
  results$p.adj <- stats::p.adjust(results$p.value, method='BH')
  # order by significance
  #results <- results[order(results$p.adj),]

  return(results)
}
I <- getSpatialPatterns(rand, weight)
range(I$observed)
head(I)
I[which(I$p.value < 0.05),]
sum(I$p.value < 0.05)/M

alpha = 0.05
lisa <- sapply(rownames(rand), function(g) {
  gexp <- rand[g,]
  Ii <- lisaTest(gexp, weight)
  lisa <- Ii$p.value
  names(lisa) <- rownames(Ii)
  sum(lisa < alpha)/length(lisa) ## percent significant
})
plot(lisa, -log10(I$p.value))
abline(h=-log10(0.05), col='red')
abline(v=0.1, col='red')
sum(I$p.value < 0.05)/M
sum(lisa > 0.1)/M
sum(I$p.value < 0.05 & lisa > 0.1)/M

