## Look at correlation between right and left side

source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

########################### Analyze old results
## pathToTable = '/n/boslfs/LABS/zhuang_lab/Everyone/180115_final_full_cluster_data/feature_table.csv'
## pathToSignal = '/n/boslfs/LABS/zhuang_lab/Everyone/180115_final_full_cluster_data/feature_signal.csv'
## pathToSignalNames = '/n/boslfs/LABS/zhuang_lab/Everyone/180115_final_full_cluster_data/signal_names.csv'
## features <- read.csv(pathToTable)
## cd <- read.csv(pathToSignal, header=FALSE)
## geneNames <- read.csv(pathToSignalNames, header=FALSE)
## colnames(cd) <- geneNames[,1]
## rownames(cd) <- features$feature_uID
## head(cd)
## save(features, cd, file='mPOA.RData')

load('mPOA.RData')

vi <- features$dataset_name=='171021_FN7_2_M22_M26'
table(vi)
features <- features[vi,]

vi <- features$sliceID == 4
table(vi)
features <- features[vi,]

vi <- features$keep_feature==1
table(vi)
features <- features[vi,]

annot <- paste0(features$cluster_name_1, features$cluster_name_2)
names(annot) <- features$feature_uID
head(annot)
table(annot)/length(annot)*100
table(annot)

pos <- features[, c('centroid_1', 'centroid_2')]
rownames(pos) <- features$feature_uID

## filter
cells.have <- features$feature_uID
annot <- annot[cells.have]
annot <- as.factor(annot)
pos <- pos[cells.have,]
cd <- cd[cells.have,]

# potential confounders
res <- features[, c("abs_volume", "backgroundScore650", "backgroundScore750", "primary_fovID")]
rownames(res) <- features$feature_uID
res <- t(res)
res <- rbind(res, libSize=colSums(cd))
res <- rbind(res, libComp=colSums(cd>0))
head(res[1:5,1:5])

plot(pos, col=annot)
plot(pos, col=map2col(scale(log10(res[1,]+1))[,1]), main='abs_volume')
plot(pos, col=map2col(scale(res[2,])[,1]), main="backgroundScore650")
plot(pos, col=map2col(scale(res[3,])[,1]), main="backgroundScore750")
plot(pos, col=res[4,], main="fovID")

mat <- normalizeCounts(t(cd))
mat <- rbind(mat, res[, colnames(mat)])
mat[1:5,1:5]
mat[is.na(mat)] <- 0

########### Run on half
plot(pos, col=annot)
range(pos$centroid_1)
median(pos$centroid_1)
pos.left <- pos[pos$centroid_1 < median(pos$centroid_1),]
pos.right <- pos[pos$centroid_1 >= median(pos$centroid_1),]

par(mfrow=c(1,2))
plot(pos.left, col=annot)
plot(pos.right, col=annot)

## run on left
adj.left <- getSpatialWeights(pos.left, k=3)
meringue.results.left <- getSpatialPatterns(mat[, rownames(pos.left)], adj.left, ncores=10)
adj.right <- getSpatialWeights(pos.right, k=3)
meringue.results.right <- getSpatialPatterns(mat[, rownames(pos.right)], adj.right, ncores=10)

## look at correlation
head(meringue.results.left)
head(meringue.results.right)

#lr <- unlist(meringue.results.left$observed)
#rr <- unlist(meringue.results.right$observed)
lp <- unlist(meringue.results.left$p.value)
rp <- unlist(meringue.results.right$p.value)

par(mfrow=c(1,1))
#plot(lr, rr)
plot(-log10(lp), -log10(rp))

########### MERingue on one cell type
table(annot)
cc <- 'OD Immature'
ct <- na.omit(names(annot)[grepl(cc,annot)])
unique(annot[ct])
plot(pos, col=c('grey', 'red')[as.factor(rownames(pos) %in% ct)],pch=16)

## run on left
ct.left <- intersect(ct, rownames(pos.left))
adj.left <- getSpatialWeights(pos.left[ct.left,], k=3)
meringue.results.left <- getSpatialPatterns(mat[, ct.left], adj.left, ncores=10)
ct.right <- intersect(ct, rownames(pos.right))
adj.right <- getSpatialWeights(pos.right[ct.right,], k=3)
meringue.results.right <- getSpatialPatterns(mat[, ct.right], adj.right, ncores=10)

## look at correlation
head(meringue.results.left)
head(meringue.results.right)

#lr <- unlist(meringue.results.left$observed)
#rr <- unlist(meringue.results.right$observed)
lp <- unlist(meringue.results.left$p.adj)
rp <- unlist(meringue.results.right$p.adj)
names(lp) <- rownames(meringue.results.left)
names(rp) <- rownames(meringue.results.right)

par(mfrow=c(1,1))
#plot(lr, rr)
plot(-log10(lp), -log10(rp))
abline(v=-log10(0.05), col='red')
abline(h=-log10(0.05), col='red')

###### ggrepel plot
pv <- data.frame(left=-log10(lp), right=-log10(rp))
sig.genes <- na.omit(rownames(pv)[pv$left > -log10(0.05) | pv$right > -log10(0.05)])
label <- rownames(pv)
label[!(label %in% sig.genes)] <- NA
library(ggplot2)
library(ggrepel)
ggplot(pv) +
  geom_point(aes(left, right), color = 'red') +
  geom_text_repel(aes(left, right, label = label)) +
  theme_classic(base_size = 16)

par(mfrow=c(1,2))
interpolate(pos[ct,], gexp=mat['Slco1a4',ct])


