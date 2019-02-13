# Try to find robust case where SpatialDE fails

source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

# Simulate Data
set.seed(1)
pos <- rbind(
 cbind(x=runif(500, 0, 100), y=runif(500, 0, 100)),
 cbind(x=runif(500, 50, 60), y=runif(500, 50, 60))
)
rownames(pos) <- paste0('cell', 1:nrow(pos))
plot(pos)

## these simulations fail for trensceek since in-homogenous poisson
#covar <- cov(t(pos))
#d <- as.matrix(dist(pos))^2
#heatmap(d)

#cd <- d[1:5,]
set.seed(1234)
cd <- matrix(runif(nrow(pos)*100), nrow=nrow(pos))*100
rownames(cd) <- rownames(pos)
cd <- t(cd)
dim(cd)
rownames(cd) <- paste0('gene', 1:nrow(cd))
plot(pos, col=map2col(cd[3,]), pch=16)

## ########## gene expression is uniformly distributed (true negative)
## set.seed(1)
## #value <- rnorm(nrow(pos))
## value <- runif(n=nrow(pos))
## names(value) <- rownames(pos)
## plot(pos, col=map2col(value), pch=16)
## v1 <- value

## m <- cov(t(pos))
## hc <- hclust(1-as.dist(m), method='ward')
## heatmap(m, ColSideColors=rainbow(2)[group], Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc))
## plot(pos, col=rainbow(2)[cutree(hc,2)])
## v2 <- m[1,]
## plot(pos, col=map2col(v2), pch=16)

## ########## gene expression is confounded with spatial pattern
## set.seed(1)
## value <- runif(nrow(pos))
## value[1:500] <- value[1:500] + rnorm(500, mean=5)
## #value <- abs(pos$x) + abs(pos$y)
## names(value) <- rownames(pos)
## plot(pos, col=map2col(value), pch=16)
## v3 <- value

## ########## gene expression
## ## set.seed(0)
## ## value <- runif(n=nrow(pos))
## ## value[pos$x > 150] <- value[pos$x > 150]*2
## ## value[pos$x < 50] <- value[pos$x < 50]*2
## ## plot(pos, col=map2col(value), pch=16)
## ## v3 <- value

## cd <- rbind(v1, v2, v3)
## head(cd)

########### MERingue
adj <- getSpatialWeights(pos, k=3)
meringue.results <- getSpatialPatterns(cd, adj)
meringue.results <- meringue.results[order(meringue.results$p.value),]
head(meringue.results)

########### SpatialDE
library(reticulate)
use_python("~/.conda/envs/testJeanFan/bin/python")
SpatialDE <- import("SpatialDE")

results = SpatialDE$run(r_to_py(as.matrix(pos)), r_to_py(as.data.frame(t(cd))))
results <- results[order(results$pval, decreasing=FALSE),]
rownames(results) <- results$g
head(results)

head(results[rownames(meringue.results),])
head(meringue.results[results$g,])

## compare
spatialde.p <- results$pval; names(spatialde.p) <- results$g
meringue.p <- meringue.results$p.val; names(meringue.p) <- rownames(meringue.results)
plot(-log10(spatialde.p), -log10(meringue.p))

g <- results$g[1]
plot(pos[,1], pos[,2], col=map2col(mat[g,]), pch=16, cex=1)

