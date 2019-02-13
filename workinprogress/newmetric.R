source('../R/helper.R')
source('../R/main.R')
source('../R/process.R')

######################################### simulate spatially distributed data
pos <- rbind(
  matrix(runif(1000, 0, 100), ncol=2),
  matrix(runif(1000, 100, 200), ncol=2)
)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- paste0('cell', 1:nrow(pos))

head(pos)

########## gene expression is uniformly distributed (true negative)
set.seed(0)
value <- runif(nrow(pos))
names(value) <- rownames(pos)
plot(pos, col=map2col(value), pch=16)

# plot
adj <- getAdj(pos, k=3)

# calculate Moran's I
library(ape)
ape::Moran.I(value, adj, alternative='greater')

# Moran's I from scratch
weight <- adj
x <- value

N <- length(x)
W <- sum(weight)
m <- mean(x)
y <- x - m
cv <- sum(weight * y %o% y)
v <- sum(y^2)
obs <- (N/W) * (cv/v)
obs

######## New correlation instead of autocorrelation metric
weight <- adj
x <- value
y <- x
plot(x,y)

N <- length(x)
W <- sum(weight)
dx <- x - mean(x)
dy <- y - mean(y)
cv <- sum(weight * dx %o% dy)
v <- sum(dx * dy)
obs <- (N/W) * (cv/v)
obs



y <- x + rnorm(length(x))/10
plot(x,y)
plot(pos, col=map2col(x), pch=16)
plot(pos, col=map2col(y), pch=16)



### tiny example
x <- c(1.5,1.5,1.5,0.1,0.1)
y <- c(1.1,1.1,0.1,1.5,1.1)
weight <- matrix(c(0,0,0,1,1,
                   0,0,0,1,1,
                   0,0,0,1,1,
                   1,1,1,0,0,
                   1,1,1,0,0), ncol=5)
rownames(weight) <- colnames(weight) <- names(x) <- names(y) <- c('A','B','C','D','E')

set.seed(0)
pos <- matrix(runif(10, 0, 100), ncol=2)
pos <- data.frame(pos)
colnames(pos) <- c('x', 'y')
rownames(pos) <- c('A','B','C','D','E')

line.col <- 'red'
col <- c('blue', 'blue', 'blue', 'green', 'green'); names(col) <- c('A','B','C','D','E')
  plot(pos, pch=16, col=col)
  idx <- which(weight>0, arr.ind = T)
  for(i in seq_len(nrow(idx))) {
    lines(
      c(pos[idx[i,1],1], pos[idx[i,2],1]),
      c(pos[idx[i,1],2], pos[idx[i,2],2]),
      col=line.col,
      lwd=weight[idx]
    )
  }

plot(pos, pch=16, col=map2col(x))
plot(pos, pch=16, col=map2col(y))

spatialCor(x, x, weight)
Moran.I(x, weight, alternative='greater')






spatialCor <- function(x, y, weight) {
    # scale weights
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight/rs

    N <- length(x)
    W <- sum(weight)
    dx <- x - mean(x)
    dy <- y - mean(y)
    cv <- sum(weight * dx %o% dy)
    v <- sqrt(sum(dx^2) * sum(dy^2))
    obs <- (N/W) * (cv/v)

    # first moment
    ei <- -1/(N - 1)

    # second moment
    W.sq <- W^2
    N.sq <- N^2
    S1 <- 0.5 * sum((weight + t(weight))^2)
    S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
    S3 <- (sum((dx*dy)^2)/N)/(v/N)^2
    S4 <- (N.sq - 3*N + 3)*S1 - N*S2 + 3*W.sq
    S5 <- (N.sq - N)*S1 - 2*N*S2 + 6*W.sq
    ei2 <- (N*S4 - S3*S5)/((N - 1)*(N - 2)*(N - 3) * W.sq)

    # standard deviation
    sdi <- sqrt(ei2 - (ei)^2)

    pv <- pnorm(obs, mean = ei, sd = sdi)
    pv <- 1 - pv

    return(list(observed = obs, expected = ei, sd = sdi, p.value = pv))
}

spatialCorPermutation <- function(x, y, weight, N=100, ncores=1) {
    # Compute Moran's I
    stat <- spatialCor(x, y, weight)$observed
    # Simulate null distribution
    sim <- unlist(parallel::mclapply(seq_len(N), function(i) {
                                spatialCor(sample(x, length(x), replace=TRUE),
                                           sample(y, length(y), replace=TRUE),
                                           weight)
                            }, mc.cores=ncores))
    sim[is.nan(sim)] <- 0
    all <- c(stat, sim)
    p.value <- mean(all >= stat, na.rm=TRUE)
    #hist(sim, sub=paste("p =", round(p.value, 4)), xlim=range(all, na.rm=TRUE))
    #abline(v = stat, col="red", lty=3, lwd=2)
    results <- unlist(data.frame('observed'=stat, 'N'=N, 'expected'=mean(all), 'sd'=sd(all), 'p.value'=p.value))
    results
}







############### Test on mPOA?
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

## limit to just MERFISH probes
vi <- c("4732456N10Rik", "Ace2", "Adora2a", "Aldh1l1", "Amigo2", "Ano3",
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
"Blank-1", "Blank-2", "Blank-3", "Blank-4", "Blank-5")
cd <- cd[, vi]

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

######## derive weight matrix between oligodendrocytes and nearest neurons
annot[grepl('Doublet', annot)] <- NA
table(annot)
## get neighbors of different cell type
cct <- na.omit(names(annot)[grepl('OD', annot)])
nct <- na.omit(names(annot)[grepl('Excitatory|Inhibitory', annot)])

k = 6
knn <- RANN::nn2(pos[nct,], pos[cct,], k=k)[[1]]
adj <- matrix(0, length(c(cct, nct)), length(c(cct, nct)))
rownames(adj) <- colnames(adj) <- c(cct, nct)
invisible(lapply(seq_len(length(cct)), function(i) {
    adj[cct[i],nct[knn[i,]]] <<- 1
}))
dim(adj)

## highly correlation due to signal bleeding though
x <- mat['Ar', c(cct, nct)]
y <- mat['Ar', c(cct, nct)]
I1 <- spatialCor(x, y, adj)
I2 <- spatialCorPermutation(x, y, adj, ncores=10)
r <- I1$obs
N <- length(c(cct, nct))
t=r/sqrt((1-r^2)/(N-2)) # correlation t statistic
pt(t, length(cct)+length(nct)-2, lower=FALSE)

rownames(mat)


# do all pairs
genes <- c('Ar', 'Egr2', 'Esr1', 'Pdgfra', 'Calcr', 'Sox4', 'Sox6', 'Sox8')
SCI <- do.call(cbind, lapply(genes, function(i) {
    print(i)
    do.call(rbind, lapply(genes, function(j) {
        print(j)
        x <- mat[i, c(cct, nct)]
        y <- mat[j, c(cct, nct)]
        spatialCor(x, y, adj)$p.value
    }))
}))
rownames(SCI) <- colnames(SCI) <- genes
SCI.melt <- reshape2::melt(SCI)
SCI.melt$value <- -log10(p.adjust(SCI.melt$value))
SCI.melt


x <- mat['Sox4', c(cct, nct)] # OD
y <- mat['Esr1', c(cct, nct)] # neuron
spatialCor(x, y, adj)

# Plot small section

foo <- pos[c(cct, nct),]
vi <- foo$centroid_1 < 1200 & foo$centroid_2 > 3000
foo <- foo[vi,]
#ccol <- rep('grey', nrow(foo)); names(ccol) <- rownames(foo)
#ccol[cct] <- 'red'
#ccol[nct] <- 'blue'
#plotNetwork(foo, weight[vi, vi], col=ccol)

plot(pos[intersect(rownames(foo), cct),], col=map2col(mat['Sox4', intersect(rownames(foo), cct)], colorRampPalette(c('green', 'white', 'orange'))(100)), pch=16)
points(pos[intersect(rownames(foo), nct),], col=map2col(mat['Esr1', intersect(rownames(foo), nct)]), pch=16)
idx <- which(weight[vi,vi]>0, arr.ind = T)
for(i in seq_len(nrow(idx))) {
    lines(
        c(foo[idx[i,1],1], foo[idx[i,2],1]),
        c(foo[idx[i,1],2], foo[idx[i,2],2]),
        col='grey',
        lwd=1
    )
}

plot(pos[intersect(rownames(foo), cct),], col=map2col(mat['Esr1', intersect(rownames(foo), cct)], colorRampPalette(c('green', 'white', 'orange'))(100)), pch=16)
points(pos[intersect(rownames(foo), nct),], col=map2col(mat['Sox4', intersect(rownames(foo), nct)]), pch=16)

plot(pos[rownames(foo),], col=map2col(mat['Esr1', rownames(foo)], colorRampPalette(c('green', 'white', 'orange'))(100)), pch=16)
plot(pos[rownames(foo),], col=map2col(mat['Sox4', rownames(foo)]), pch=16)
plot(pos, col=map2col(mat['Esr1',], colorRampPalette(c('green', 'white', 'orange'))(100)), pch=16)
plot(pos, col=map2col(mat['Sox4',]), pch=16)








## # Cross correlation: https://arxiv.org/pdf/1503.02908.pdf
## spatialXCor <- function(x, y, weight) {
##     # standarization
##     dx <- (x - mean(x))/sd(x)
##     dy <- (y - mean(y))/sd(y)

##     # unitization
##     rs <- rowSums(weight)
##     rs[rs == 0] <- 1
##     W <- weight/rs

##     # measurement for spatial cross correlation
##     SCI <- t(dx) %*% W %*% dy
##     return(SCI)
## }
## spatialXCor(x, x, weight)
## Moran.I(x, weight)





########## Use simulation
library('trendsceek')

##create synthetic dataset
#pp = sim_pois(300)
set.seed(0)
win_len = 1
lambda_int = function(x,y) {1000 * exp(-3*y)} # inhomogenous poisson
#lambda_int = 300
pp = spatstat::rpoispp(lambda_int, win = spatstat::owin(c(0, win_len), c(0, win_len)))
low_expr = c(10, 10)
high_expr = c(20, 50)
pp = add_markdist_hotspot(pp, low_expr, high_expr, hotspot_size=0.25)
pp = add_markdist_streak(pp, low_expr, high_expr, streak_height_frac=0.25)
pp = add_markdist_step(pp, low_expr, high_expr, step_border=0.25)
pos <- cbind(pp$x, pp$y)
cd <- t(pp$marks) + rnorm(nrow(pos))*5
#cd <- t(pp$marks)
dim(cd)
rownames(pos) <- colnames(cd)
plot(pos, col=map2col(cd[1,]), pch=16)
plot(pos, col=map2col(cd[3,]), pch=16)
plot(pos, col=map2col(cd[5,]), pch=16)

adj <- getSpatialWeights(pos, k=3)

Moran.I(cd[1,], adj, alternative='greater')
spatialCor(cd[1,], cd[1,], adj)

moranPermutationTest(cd[1,], adj, N=100)
spatialCorPermutation(cd[1,], cd[1,], adj, N=1000)

cor(cd[1,], cd[2,])
spatialCor(cd[1,], cd[2,], adj)
spatialCorPermutation(cd[1,], cd[2,], adj, N=1000)

cor(cd[1,], cd[6,])
spatialCor(cd[1,], cd[6,], adj)
spatialCorPermutation(cd[1,], cd[6,], adj, N=1000)


################################### simulate perfect data
x <- rnorm(1000, mean=10, sd=1)
y <- x + rnorm(1000)
plot(x,y) # highly correlated

pos <- cbind(x, y)
rownames(pos) <- paste0('cell', 1:nrow(pos))
head(pos)

par(mfrow=c(1,2))
y <- rnorm(1000, mean=10, sd=1)
plot(pos, col=map2col(x), pch=16)
plot(pos, col=map2col(y), pch=16)

adj <- getSpatialWeights(pos, k=3)
#plotNetwork(pos, adj)

unlist(spatialCor(x, y, adj))
spatialCorPermutation(x, y, adj, N=1000)

unlist(spatialCor(x, x, adj))
spatialCorPermutation(x, x, adj, N=1000)
unlist(Moran.I(x, adj, alternative='greater'))

unlist(spatialCor(y, y, adj))
spatialCorPermutation(y, y, adj, N=1000)
unlist(Moran.I(y, adj, alternative='greater'))
