source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

########################### Analyze old results
features <- read.csv('/n/boslfs/LABS/zhuang_lab/Everyone/180115_final_full_cluster_data/feature_table.csv')
vi <- features$dataset_name=='171021_FN7_2_M22_M26'
table(vi)
features <- features[vi,]

vi <- features$sliceID == 5
table(vi)
features <- features[vi,]

annot <- paste0(features$cluster_name_1, features$cluster_name_2)
names(annot) <- features$feature_uID
head(annot)
table(annot)/length(annot)*100
table(annot)

pos <- features[, c('centroid_1', 'centroid_2')]
rownames(pos) <- features$feature_uID

## read codebook
codebook <- read.csv('../../../mPOA/M22E1_codebook.csv', skip=3)
## load old counts matrix
library(Matrix)
cellCounts1.jeff <- Matrix(as.matrix(read.csv("/n/boslfs/LABS/zhuang_lab/User/jeffmoffitt/normalized_data/171021_FN7_2_M22_M26/reports/countsPerCellExactIn.csv", header=FALSE)), sparse=TRUE)
cellCounts2.jeff <- Matrix(as.matrix(read.csv("/n/boslfs/LABS/zhuang_lab/User/jeffmoffitt/normalized_data/171021_FN7_2_M22_M26/reports/countsPerCellCorrectedIn.csv", header=FALSE)), sparse=TRUE)
cellCounts.jeff = cellCounts1.jeff + cellCounts2.jeff
dim(cellCounts.jeff)
rm(cellCounts1.jeff)
rm(cellCounts2.jeff)
gc()
rownames(cellCounts.jeff) <- codebook$name

cellNames <- read.csv("/n/boslfs/LABS/zhuang_lab/User/jeffmoffitt/normalized_data/171021_FN7_2_M22_M26/reports/featureNames.csv", header=FALSE)
colnames(cellCounts.jeff) <- cellNames[,1]

## filter
vi <- annot != ""
annot <- annot[vi]
cells.have <- intersect(names(annot), rownames(pos))

annot <- annot[cells.have]
annot <- as.factor(annot)
pos <- pos[cells.have,]
cd <- cellCounts.jeff[, cells.have]

plot(pos, col=annot)

## limit to subtype
levels(annot)
pos.sub <- pos[annot=="Excitatory" ,]
cd.sub <- cd[, annot=="Excitatory" ]

########### MERingue
adj <- getSpatialWeights(pos.sub, k=6)
mat <- normalizeCounts(cd.sub)
meringue.results <- getSpatialPatterns(mat, adj, ncores=4)
head(meringue.results)

########### SpatialDE
library(reticulate)
use_python("~/.conda/envs/testJeanFan/bin/python")
SpatialDE <- import("SpatialDE")

results = SpatialDE$run(r_to_py(as.matrix(pos.sub)), r_to_py(as.data.frame(t(as.matrix(mat)))))
results <- results[order(results$qval, decreasing=FALSE),]
head(results)

spatialde.p <- results$qval; names(spatialde.p) <- results$g
meringue.p <- meringue.results$p.adj; names(meringue.p) <- rownames(meringue.results)
plot(-log10(spatialde.p), -log10(meringue.p), col=rainbow(2)[as.factor(grepl('Blank', names(spatialde.p)))], pch=16)
abline(v=-log10(0.05))
abline(h=-log10(0.05))

rownames(meringue.results)[meringue.results$p.adj < 0.05]
results$g[results$qval < 0.05]
g <- 'Irs4'
plot(pos.sub[,1], pos.sub[,2], col=map2col(mat[g,]), pch=16, cex=1)

########## Trendsceeek
library('trendsceek')
pp = pos2pp(pos)
log.fcn = log10
pp = set_marks(pp, as.matrix(cd), log.fcn = log.fcn)
pp2plot = pp_select(pp)
##set parameters
nrand = 100
ncores = 1
##run
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

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
