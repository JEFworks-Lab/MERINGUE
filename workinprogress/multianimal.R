## Look at correlation across animals
source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

########################### Load data

load('mPOA.RData')

getData <- function(features, cd, s1) {
    vi <- features$dataset_name==s1
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
    #plot(pos, col=map2col(scale(log10(res[1,]+1))[,1]), main='abs_volume')
    #plot(pos, col=map2col(scale(res[2,])[,1]), main="backgroundScore650")
    #plot(pos, col=map2col(scale(res[3,])[,1]), main="backgroundScore750")
    #plot(pos, col=res[4,], main="fovID")

    mat <- normalizeCounts(t(cd))
    #mat <- rbind(mat, res[, colnames(mat)])
    mat[1:5,1:5]
    mat[is.na(mat)] <- 0

    return(list(
        'mat'=mat,
        'res'=res,
        'annot'=annot,
        'pos'=pos))
}


unique(features$dataset_name)
s1 <- getData(features, cd, '171021_FN7_2_M22_M26')
s2 <- getData(features, cd, '170921_FN4_2_M22_M26')

########### Run on both mice
plot(s1$pos, col=s1$annot)
plot(s2$pos, col=s2$annot)

adj.s1 <- getSpatialWeights(s1$pos, k=3)
meringue.results.s1 <- getSpatialPatterns(s1$mat, adj.s1, ncores=20)
adj.s2 <- getSpatialWeights(s2$pos, k=3)
meringue.results.s2 <- getSpatialPatterns(s2$mat, adj.s2, ncores=20)

lp <- unlist(meringue.results.s1$p.adj)
rp <- unlist(meringue.results.s2$p.adj)
names(lp) <- rownames(meringue.results.s1)
names(rp) <- rownames(meringue.results.s2)

## plot
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

########### one cell type
table(s1$annot)
cc <- 'OD Mature 2'

ct1 <- na.omit(names(s1$annot)[grepl(cc, s1$annot)])
plot(s1$pos, col=c('grey', 'red')[as.factor(rownames(s1$pos) %in% ct1)],pch=16)

ct2 <- na.omit(names(s2$annot)[grepl(cc, s2$annot)])
plot(s2$pos, col=c('grey', 'red')[as.factor(rownames(s2$pos) %in% ct2)],pch=16)

adj.s1 <- getSpatialWeights(s1$pos[ct1,], k=3)
meringue.results.s1 <- getSpatialPatterns(s1$mat[, ct1], adj.s1, ncores=20)
adj.s2 <- getSpatialWeights(s2$pos[ct2,], k=3)
meringue.results.s2 <- getSpatialPatterns(s2$mat[, ct2], adj.s2, ncores=20)

lp <- unlist(meringue.results.s1$p.adj)
rp <- unlist(meringue.results.s2$p.adj)
names(lp) <- rownames(meringue.results.s1)
names(rp) <- rownames(meringue.results.s2)

## plot
pv <- data.frame(s1=-log10(lp), s2=-log10(rp))
sig.genes <- na.omit(rownames(pv)[pv$left > -log10(0.05) | pv$right > -log10(0.05)])
label <- rownames(pv)
label[!(label %in% sig.genes)] <- NA
library(ggplot2)
library(ggrepel)
ggplot(pv) +
  geom_point(aes(left, right), color = 'red') +
  geom_text_repel(aes(left, right, label = label)) +
  theme_classic(base_size = 16)

g <- 'Opalin'
plot(s1$pos[ct1,], col=map2col(s1$mat[g, ct1]), pch=16)
plot(s2$pos[ct2,], col=map2col(s2$mat[g, ct2]), pch=16)
