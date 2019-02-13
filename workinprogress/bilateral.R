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

vi <- features$sliceID == 5
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

cd <- cd[, c("4732456N10Rik", "Ace2", "Adora2a", "Aldh1l1", "Amigo2", "Ano3",
"Aqp4", "Ar", "Arhgap36", "Avpr1a", "Avpr2", "Baiap2", "Bdnf",
"Bmp7", "Brs3", "Calcr", "Cbln1", "Cbln2", "Cckar", "Cckbr",
"Ccnd2", "Cd24a", "Cdkn1a", "Cenpe", "Chat", "Coch", "Col25a1",
"Cplx3", "Cpne5", "Creb3l1", "Crhbp", "Crhr1", "Crhr2", "Cspg5",
"Cxcl14", "Cyp19a1", "Cyp26a1", "Cyr61", "Dgkk", "Ebf3", "Egr2",
"Ermn", "Esr1", "Etv1", "Fbxw13", "Fezf1", "Fn1", "Fst", "Gabra1",
"Gabrg1", "Gad1", "Galr1", "Galr2", "Gbx2", "Gda", "Gem", "Gjc3",
"Glra3", "Gpr165", "Greb1", "Grpr", "Htr2c", "Igf1r", "Igf2r",
"Irs4", "Isl1", "Kiss1r", "Klf4", "Lepr", "Lmod1", "Lpar1", "Man1a",
"Mc4r", "Mki67", "Mlc1", "Myh11", "Ndnf", "Ndrg1", "Necab1",
"Nos1", "Npas1", "Npy1r", "Npy2r", "Ntng1", "Ntsr1", "Nup62cl",
"Omp", "Onecut2", "Opalin", "Oprd1", "Oprk1", "Oprl1", "Oxtr",
"Pak3", "Pcdh11x", "Pdgfra", "Pgr", "Plin3", "Pnoc", "Pou3f2",
"Prlr", "Ramp3", "Rgs2", "Rgs5", "Rnd3", "Rxfp1", "Scgn", "Selplg",
"Sema3c", "Sema4d", "Serpinb1b", "Serpine1", "Sgk1", "Slc15a3",
"Slc17a6", "Slc17a7", "Slc17a8", "Slc18a2", "Slco1a4", "Sox4",
"Sox6", "Sox8", "Sp9", "Synpr", "Syt2", "Syt4", "Sytl4", "Tacr1",
"Tacr3", "Tiparp", "Tmem108", "Traf4", "Trhr", "Ttn", "Ttyh2",
"Blank-1", "Blank-2", "Blank-3", "Blank-4", "Blank-5")]

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
#head(meringue.results.left)
#head(meringue.results.right)

#lr <- unlist(meringue.results.left$observed)
#rr <- unlist(meringue.results.right$observed)
lp <- unlist(meringue.results.left$p.adj)
rp <- unlist(meringue.results.right$p.adj)

par(mfrow=c(1,1))
#plot(lr, rr)
plot(-log10(lp), -log10(rp))
abline(v=-log10(0.05), col='red')
abline(h=-log10(0.05), col='red')

########### MERingue on one cell type
table(annot)
cc <- 'Inhibitory'
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
#head(meringue.results.left)
#head(meringue.results.right)

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
g <- 'Egr2'
interpolate(pos[ct,], gexp=mat[g,ct])


