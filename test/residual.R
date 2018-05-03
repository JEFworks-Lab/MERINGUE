## residual testing on mPOA data

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

vi <- features$keep_feature==1
features <- features[vi,]

head(features)

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
#cells.have <- intersect(names(annot), rownames(pos))
cells.have <- features$feature_uID

annot <- annot[cells.have]
annot <- as.factor(annot)
pos <- pos[cells.have,]
cd <- cellCounts.jeff[, cells.have]

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

########### MERingue on confounders
ct <- names(annot)[annot=='OD Mature 2']
adj <- getSpatialWeights(pos[ct,], k=3)
mat <- normalizeCounts(cd[, ct])
mat <- rbind(mat, res[, colnames(mat)])
#mat <- res
meringue.results <- getSpatialPatterns(mat, adj, ncores=1)
head(meringue.results)
meringue.results[rownames(res),]

########### look at residuals
#So, when a problem has a spatial component, we should:
#Run the non-spatial regression
#Test the regression residuals for spatial autocorrelation, using Moran's I or some other index
#If no significant spatial autocorrelation exists, STOP. Otherwise, if the spatial dependencies are significant, use a special model which takes spatial dependencies into account.
table(annot)
m <- t(as.matrix(mat[, ct]))
#m[m==0] <- NA
hh <- as.data.frame(m)
class(hh)
f1 <- abs_volume ~ backgroundScore650
plot(m[,'backgroundScore650'], m[,'abs_volume'])
m1 <- lm(f1, data=hh, na.action=na.omit)
summary(m1)
#plot(m1)
moranTest(x=m[names(m1$residuals),'abs_volume'], w=adj[names(m1$residuals),names(m1$residuals)]) # spatial dependency exists...
moranTest(x=m1$residuals, w=adj[names(m1$residuals),names(m1$residuals)]) # spatial dependency exists...
#plot(pos[names(m1$residuals),], col=map2col(t(mat)[names(m1$residuals),'Ar']))
plot(pos[names(m1$residuals),], col=map2col(m[names(m1$residuals),'abs_volume']))
plot(pos[names(m1$residuals),], col=map2col(m1$residuals))

plot(pos[names(m1$residuals),], col=map2col(t(mat)[names(m1$residuals),'backgroundScore750']))
moranTest(x=t(mat)[names(m1$residuals),'backgroundScore750'], w=adj[names(m1$residuals),names(m1$residuals)]) # spatial dependency exists...

plot(pos[names(m1$residuals),], col=map2col(t(mat)[names(m1$residuals),'Blank-5']))
moranTest(x=t(mat)[names(m1$residuals),'Blank-5'], w=adj[names(m1$residuals),names(m1$residuals)]) # spatial dependency exists...

plot(pos, col=map2col(t(mat)[,'Blank-1']), pch=16)
moranTest(x=t(mat)[,'Blank-1'], w=adj) # spatial dependency exists...

########### See if subtype patterns can be explained by general patterns in Ar expression
interpolate <- function(pos, gexp, binSize=100, cex=1, col=colorRampPalette(c('blue', 'white', 'red'))(100), plot=TRUE, ...) {
  z <- gexp
  x <- pos[,1]
  y <- pos[,2]
  int <- akima::interp(x, y, z, nx=binSize, ny=binSize, linear=FALSE)
  if(plot) {
    plot(pos, col=map2col(z), pch=16, cex=cex, axes=FALSE, frame.plot=TRUE, xlab=NA, ylab=NA, ...)
    image(int, col=col, axes=FALSE, frame.plot=TRUE, ...)
  }
  return(int)
}

mat <- normalizeCounts(cd)
g <- 'Slc18a2'

## get non-oligodendrocyte neighbors
annot[grepl('Doublet', annot)] <- NA
ct <- na.omit(names(annot)[annot=='Inhibitory'])
nct <- na.omit(names(annot)[annot!='Inhibitory'])
#nct <- na.omit(names(annot)[annot=='OD Mature 2'])
k = 3
knn <- RANN::nn2(pos[nct,], pos[ct,], k=k)[[1]]
adj <- matrix(0, nrow(pos), nrow(pos))
rownames(adj) <- colnames(adj) <- rownames(pos)
invisible(lapply(seq_len(length(ct)), function(i) {
    adj[ct[i],nct[knn[i,]]] <<- 1
}))

par(mfrow=c(2,2), mar=rep(1,4))
#interpolate(pos, mat[g,], cex=0.5, binSize=100)
interpolate(pos[ct,], mat[g,ct], cex=0.5, binSize=100)
interpolate(pos[nct,], mat[g,nct], cex=0.5, binSize=100)
moranTest(x=mat['Slc18a2',nct], w=getAdj(pos[nct,], k=3))

# plot expression
g <- 'Slc18a2'
x <- mat[g, ct]
y <- do.call(rbind, lapply(seq_len(length(ct)), function(i) mat[g,nct[knn[i,]]]))
plot(x,y[,3])

plotNetwork <- function(pos, adj, col='black', line.col='red', line.power=1, ...) {
  plot(pos, pch=16, col=col)
  idx <- which(adj>0, arr.ind = T)
  for(i in seq_len(nrow(idx))) {
    lines(
      c(pos[idx[i,1],1], pos[idx[i,2],1]),
      c(pos[idx[i,1],2], pos[idx[i,2],2]),
      col=line.col,
      lwd=adj[idx]^line.power
    )
  }
}

col <- c('black', 'green')[as.factor(rownames(pos) %in% ct)]
names(col) <- rownames(pos)
plotNetwork(pos[1:100,], adj[1:100,1:100], col=col)

## make adjacency matrix square
moranTest(x=mat['Slc18a2',], w=adj)

