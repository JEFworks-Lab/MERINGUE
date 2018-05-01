source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

# Data
cd <- read.csv("../../../MERingue/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv", row.names=1)
pos <- read.csv("../../../MERingue/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv", row.names=1)
cd <- cd[rownames(pos),]
sample_info <- pos
pos <- pos[,1:2]
counts <- cleanCounts(t(as.matrix(cd)), min.reads=100, min.detected=10)
mat <- normalizeCounts(counts)

# test on subset of genes
genes <- rownames(mat)[1:1000]

########### MERingue
adj <- getSpatialWeights(pos, k=6)
meringue.results <- getSpatialPatterns(mat[genes,], adj)
head(meringue.results)

########### SpatialDE
library(reticulate)
use_python("~/.conda/envs/testJeanFan/bin/python")

SpatialDE <- import("SpatialDE")
NaiveDE <- import("NaiveDE")

norm_expr = t(NaiveDE$stabilize(r_to_py(t(as.data.frame(cd)))))
resid_expr = t(NaiveDE$regress_out(r_to_py(sample_info), r_to_py(t(norm_expr)), 'np.log(total_counts)'))
rownames(resid_expr) <- rownames(cd)
colnames(resid_expr) <- colnames(cd)

results = SpatialDE$run(r_to_py(as.matrix(pos)), r_to_py(as.data.frame(resid_expr[, genes])))

results <- results[order(results$qval, decreasing=FALSE),]
head(results)

########## Trendsceeek
library('trendsceek')

pp = pos2pp(pos)
log.fcn = log10
pp = set_marks(pp, as.matrix(counts[genes,]), log.fcn = log.fcn)

pp2plot = pp_select(pp)

##set parameters
nrand = 100
ncores = 1
##run
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

head(trendstat_list[['supstats_wide']])
extract_sig_genes

supstats_list = trendstat_list[["supstats"]]
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

######### Restrict to 'cell type'

d <- cor(as.matrix(mat))
hc <- hclust(1-as.dist(d), method='ward.D')
plot(hc)
groups <- cutree(hc, 2)
plot(pos[,1], pos[,2], col=map2col(groups), pch=16, cex=1)

sub <- names(groups)[groups==2]

## spatialde
results.sub = SpatialDE$run(r_to_py(as.matrix(pos[sub,])), r_to_py(as.data.frame(resid_expr[sub, genes])))
results.sub <- results.sub[order(results.sub$qval, decreasing=FALSE),]
head(results.sub)

adj <- getSpatialWeights(pos[sub,], k=3)
meringue.results.sub <- getSpatialPatterns(mat[genes,sub], adj)
head(meringue.results.sub)

spatialde.p <- results.sub$qval; names(spatialde.p) <- results.sub$g
meringue.p <- meringue.results.sub$p.adj; names(meringue.p) <- rownames(meringue.results.sub)
plot(-log10(spatialde.p), -log10(meringue.p))

g <- names(meringue.p)[4]
plot(pos[sub,1], pos[sub,2], col=map2col(mat[g,sub]), pch=16, cex=1)
