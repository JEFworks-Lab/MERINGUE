## residual testing on mPOA data
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

########### MERingue on one cell type
table(annot)
cc <- 'Astrocyte'
ct <- na.omit(names(annot)[grepl(cc,annot)])
unique(annot[ct])
plot(pos, col=c('grey', 'red')[as.factor(rownames(pos) %in% ct)],pch=16)
adj <- getSpatialWeights(pos[ct,], k=3)
meringue.results <- getSpatialPatterns(mat[, ct], adj, ncores=1)
head(meringue.results, n=20)
meringue.results[rownames(res),]

########## Look for potential confounders
library(parallel)
foo <- mclapply(rownames(mat), function(g) {
    print(g)

    #plot(pos, col=map2col(scale(mat[g,])[,1]), pch=16)
    #plot(pos[ct,], col=map2col(scale(mat[g,ct])[,1]), pch=16)

    ## get neighbors of different cell type
    cct <- na.omit(names(annot)[grepl(cc, annot)])
    nct <- na.omit(names(annot)[!grepl(cc, annot)])

    par(mfrow=c(2,2), mar=rep(1,4))
    #interpolate(pos, mat[g,], cex=0.5, binSize=100)
    #interpolate(pos[cct,], mat[g,cct], cex=0.5, binSize=100)
    #interpolate(pos[nct,], mat[g,nct], cex=0.5, binSize=100)
    r1 <- moranTest(x=mat[g,cct], w=getAdj(pos[cct,], k=3))
    r2 <- moranTest(x=mat[g,nct], w=getAdj(pos[nct,], k=3))

    k = 3
    knn <- RANN::nn2(pos[nct,], pos[cct,], k=k)[[1]]
    adj <- matrix(0, nrow(pos), nrow(pos))
    rownames(adj) <- colnames(adj) <- rownames(pos)
    invisible(lapply(seq_len(length(cct)), function(i) {
        adj[cct[i],nct[knn[i,]]] <<- 1
    }))
    # can expression pattern be explained by neighboring non-excitatory cell types?
    r3 <- moranTest(x=mat[g,], w=adj)

    cbind(r1,r2,r3)
}, mc.cores=10)
names(foo) <- rownames(mat)
foo.p <- do.call(rbind, lapply(foo, function(x) x['p.value',]))
head(foo.p)

plot(abs(qnorm(foo.p[,1], lower.tail=FALSE)), abs(qnorm(foo.p[,3], lower.tail=FALSE)))

# can be explained by neighbors
vi <- names(which(foo.p[,1]<0.05 & foo.p[,3]<0.05))
vi1 <- vi

## cannot be explained by neighbors
vi2 <- names(which(foo.p[,1]<0.05 & foo.p[,3]>0.05))

vi <- vi2
vi <- intersect(vi, na.omit(rownames(meringue.results)[meringue.results$p.adj<0.05]))
foo.p[vi,]
meringue.results[vi,]

#p.new <- unlist(lapply(foo, function(bar) { pnorm(q=bar['observed', 'r1'], bar['expected', 'r3'], bar['sd', 'r3'], lower.tail=FALSE) }))
#p.new <- stats::p.adjust(p.new)
#vi3 <- names(which(p.new < 0.05))

g <- 'Ebf3'
foo[[g]]
#bar <- foo[[g]]
meringue.results[g,]

cct <- na.omit(names(annot)[grepl(cc, annot)])
nct <- na.omit(names(annot)[!grepl(cc, annot)])
par(mfrow=c(2,2), mar=rep(1,4))
interpolate(pos[cct,], mat[g,cct], cex=0.5, binSize=100)
interpolate(pos[nct,], mat[g,nct], cex=0.5, binSize=100)

# plot expression
## x <- mat[g, cct]
## y <- do.call(rbind, lapply(seq_len(length(cct)), function(i) mat[g,nct[knn[i,]]]))
## plot(x,rowMeans(y))
## col <- c('black', 'green')[as.factor(rownames(pos) %in% ct)]
## names(col) <- rownames(pos)
## plotNetwork(pos[1:100,], adj[1:100,1:100], col=col)


