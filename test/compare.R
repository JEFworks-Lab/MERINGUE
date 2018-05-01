source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

# Simulate Data
set.seed(1)
#pos <- rbind(
#  cbind(x=rnorm(500, 0, 10), y=rnorm(50, 0, 10)),
#  cbind(x=rnorm(50, 100, 10), y=rnorm(50, 100, 10))
#)+1000
## pos <- rbind(
##   cbind(x=rnorm(500, 0, 100), y=1:500),
##   cbind(x=rnorm(50, 0, 100), y=(1:50)+500)
## )+1000
#pos <- rbind(
#  matrix(runif(100, 0, 100), ncol=2),
#  matrix(runif(100, 100, 200), ncol=2)
#)
#pos[,1] <- sort(pos[,1])
#pos[,2] <- sort(pos[,2])
pos <- rbind(
    cbind(rep(1:10, 50), 1:500),
    cbind(rep(1:10, 50), 1:100)
)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- paste0('cell', 1:nrow(pos))
head(pos)
group <- c(rep('g1', 500), rep('g2', 50))
names(group) <- rownames(pos)
group <- factor(group)
plot(pos, col=rainbow(2)[group])

########## gene expression is uniformly distributed (true negative)
set.seed(1)
#value <- rnorm(nrow(pos))
value <- runif(n=nrow(pos))
names(value) <- rownames(pos)
plot(pos, col=map2col(value), pch=16)
v1 <- value

m <- cov(t(pos))
hc <- hclust(1-as.dist(m), method='ward')
heatmap(m, ColSideColors=rainbow(2)[group], Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc))
plot(pos, col=rainbow(2)[cutree(hc,2)])
v2 <- m[1,]
plot(pos, col=map2col(v2), pch=16)

########## gene expression is confounded with spatial pattern
set.seed(1)
value <- runif(nrow(pos))
value[1:500] <- value[1:500] + rnorm(500, mean=5)
#value <- abs(pos$x) + abs(pos$y)
names(value) <- rownames(pos)
plot(pos, col=map2col(value), pch=16)
v3 <- value

########## gene expression
## set.seed(0)
## value <- runif(n=nrow(pos))
## value[pos$x > 150] <- value[pos$x > 150]*2
## value[pos$x < 50] <- value[pos$x < 50]*2
## plot(pos, col=map2col(value), pch=16)
## v3 <- value

cd <- rbind(v1, v2, v3)
head(cd)

########### MERingue
adj <- getSpatialWeights(pos, k=3)
meringue.results <- getSpatialPatterns(cd, adj)
head(meringue.results)

########### SpatialDE
library(reticulate)
use_python("~/.conda/envs/testJeanFan/bin/python")
SpatialDE <- import("SpatialDE")

results = SpatialDE$run(r_to_py(as.matrix(pos)), r_to_py(as.data.frame(t(cd))))
results <- results[order(results$qval, decreasing=FALSE),]
head(results)

########## Trendsceeek
library('trendsceek')
pp = pos2pp(pos)
log.fcn = log10
pp = set_marks(pp, 10^as.matrix(cd), log.fcn = log.fcn)
pp2plot = pp_select(pp)
##set parameters
nrand = 100
ncores = 1
##run
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

supstats_list = trendstat_list[["supstats"]]
extract_sig_genes(trendstat_list, alpha = 0.1)
sig.list = do.call(cbind, lapply(supstats_list, function(j.df) {
    x <- j.df[, "min.pval"]
    names(x) <- j.df[, 'gene']
    x
}))
head(sig.list)

########## Compare

spatialde.p <- results$qval; names(spatialde.p) <- results$g
meringue.p <- meringue.results$p.adj; names(meringue.p) <- rownames(meringue.results)
trendsceek.p <- sig.list[, 'Emark']
plot(-log10(spatialde.p), -log10(meringue.p))
plot(-log10(spatialde.p), -log10(trendsceek.p))
plot(-log10(meringue.p), -log10(trendsceek.p))
plot(spatialde.p, meringue.p)
plot(spatialde.p, trendsceek.p)
plot(meringue.p, trendsceek.p)
cor.test(spatialde.p, meringue.p)

g <- results$g[1]
plot(pos[,1], pos[,2], col=map2col(mat[g,]), pch=16, cex=1)
plot(pos[,1], pos[,2], col=map2col(resid_expr[,g]), pch=16, cex=1)

