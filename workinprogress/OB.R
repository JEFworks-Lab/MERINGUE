# Sample analysis using OB data
set.seed(0)
library(MERingue)
data(mOB)

# Clean and normalize data
pos <- mOB$pos
cd <- mOB$counts
dim(cd)
par(mfrow=c(2,1), mar=rep(2,4))
hist(log10(Matrix::colSums(cd)+1), breaks=20)
abline(v = log10(10+1), col='red', lty=2)
hist(log10(Matrix::rowSums(cd)+1), breaks=20)
abline(v = log10(10+1), col='red', lty=2)

counts <- cleanCounts(cd, min.reads=100, min.lib.size=100)
dim(counts)
par(mfrow=c(2,1), mar=rep(2,4))
hist(log10(Matrix::colSums(counts)+1), breaks=20)
hist(log10(Matrix::rowSums(counts)+1), breaks=20)

pos <- pos[colnames(counts),]
mat <- normalizeCounts(counts, log=FALSE)

# Get transcriptional clusters
pcs.info <- prcomp(t(log10(as.matrix(mat)+1)), center=TRUE)
names(pcs.info)
par(mfrow=c(1,1), mar=rep(5,4))
plot(pcs.info$sdev[1:10], type='l', main='PCs')
nPcs <- 5
abline(v=nPcs, col='red', lty=2)
pcs <- pcs.info$x[,1:nPcs]
head(pcs)
d <- dist(pcs, method='euc')
emb <- Rtsne::Rtsne(d,
                    is_distance=TRUE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)

library(igraph)
library(RANN)
k = 30
nn = nn2(as.matrix(pcs), k = k)
nn.df = data.frame(from = rep(1:nrow(nn$nn.idx), k),
                   to = as.vector(nn$nn.idx),
                   weight = 1/(1 + as.vector(nn$nn.dists)))
nw.norm = graph_from_data_frame(nn.df, directed = FALSE)
nw.norm = simplify(nw.norm)
lc.norm = cluster_louvain(nw.norm)
com = as.factor(membership(lc.norm))
names(com) <- rownames(pcs)
head(com)

par(mfrow=c(2,1), mar=rep(2,4))
plotEmbedding(emb, groups=com, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='topleft')
plotEmbedding(pos, groups=com, cex=1, xlab=NA, ylab=NA)

annot <- as.character(com); names(annot) <- names(com)
annot[com==4] <- '1: Granular Cell Layer'
annot[com==1] <- '2: Mitral Cell Layer'
annot[com==3] <- '3: Outer Plexiform Layer'
annot[com==2] <- '4: Glomerular Layer'
annot[com==5] <- '5: Olfactory Nerve Layer'
annot <- as.factor(annot)

par(mfrow=c(2,1), mar=rep(5,4))
plotEmbedding(emb, groups=annot, show.legend=TRUE, xlab=NA, ylab=NA, legend.x='topleft')
plotEmbedding(pos, groups=annot, cex=1, xlab=NA, ylab=NA)

#################### Identify differentially expressed genes
dg <- getDifferentialGenes(as.matrix(mat), annot)
dg.sig <- lapply(dg, function(x) {
  x <- x[x$p.adj < 0.05,]
  #x <- x[p.adjust(x$p.value, method="BH") < 0.05,]
  x <- na.omit(x)
  x <- x[x$highest,]
  rownames(x)
})
print(lapply(dg.sig, length))

dg.genes <- unlist(dg.sig)
length(dg.genes)
ggroup <- unlist(sapply(names(dg.genes), function(x) strsplit(x, 'Layer')[[1]][1]))
names(ggroup) <- dg.genes
ggroup <- factor(ggroup)

library(RColorBrewer)
ccol <- rainbow(length(levels(annot)))[annot]
names(ccol) <- names(annot)
gcol <- rainbow(length(levels(ggroup)), v=0.5)[ggroup]
names(gcol) <- names(ggroup)
m <- as.matrix(mat[dg.genes,names(sort(annot))])
m <- t(scale(t(m)))
range(m)
m[m < -2.5] <- -2.5
m[m > 2.5] <- 2.5

gcol <- gcol[rownames(m)]
ccol <- ccol[colnames(m)]

# pick out a few genes to highlight in heatmap
marker.genes <- c("Penk", "Doc2g", "Kctd12",
                  "Krt17", "Gas6", "Sparc", "Vim",
                  "Mucl1", "Pip", "Fn1", "Igfbp5",
                  "Scgb2a2", "Peg10", "Areg", "Mmp14", "Dcn",
                  "Camk4", "Th", "Vip",
                  "Slc17a7", "Reln", "Cdhr1", "Sv2b", "Shisa3", "Plcxd2", "Nmb", "Uchl1", "Rcan2")
rownames(m)[!(rownames(m) %in% marker.genes)] <- NA

library(gplots)
heatmap.2(m, trace="none", Colv=NA, Rowv=NA, labCol=NA,
          ColSideColors=ccol,
          RowSideColors=gcol,
          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)


#################### Identify spatially clustered genes
#voronoiAdjacency
par(mfrow=c(1,1))
w <- voronoiAdjacency(pos, filterDist=2.5, plot=TRUE)
dim(w)

#getSpatialPatterns
m <- mat
dim(m)
I <- getSpatialPatterns(m, w)
range(I$observed)
head(I)

#getMinPercentCells
getMinPercentCells <- function(weight, mat, alpha=0.05, adjustPv = TRUE, M=10, plot=TRUE) {

  lisa <- unlist(lapply(seq_len(M), function(i) {
    ## show shuffle
    set.seed(i)
    rand <- mat[,sample(ncol(mat))]
    colnames(rand) <- colnames(mat)

    ## assess significance of randomly permuted genes
    I <- getSpatialPatterns(rand, weight)
    if(adjustPv) {
      pv <- I$p.adj
    } else {
      pv <- I$p.value
    }
    falsePositives <- rownames(I)[pv < alpha]
    print(falsePositives)

    ## get lisa
    lisa <- sapply(rownames(rand), function(g) {
      gexp <- rand[g,]
      Ii <- lisaTest(gexp, weight)
      lisa <- Ii$p.value
      names(lisa) <- rownames(Ii)
      sum(lisa < alpha)/length(lisa) ## percent significant
    })
    hist(lisa, breaks=20)
    abline(v=0.05, col='red')

    par(mfrow=c(5,2))
    g <- falsePositives[1]
    gexp <- scale(mat[g,])[,1]
    gexp[gexp > 2] <- 2
    gexp[gexp < -2] <- -2
    #invisible(interpolate(pos[colnames(m),], m[g,], main=g, zlim=c(-2,2)))
    plot(pos, col=MERingue:::map2col(gexp), pch=16, cex=0.75, axes=FALSE, frame.plot=TRUE, main=g)
    slisa <- MERingue:::getSignedLisa(gexp, w, plot=FALSE)
    #slisa <- scale(slisa)[,1]
    # slisa[slisa < -1] <- -1
    #slisa[slisa > 1] <- 1
    plot(pos, col=MERingue:::map2col(slisa, pal=colorRampPalette(c('darkgreen', 'white', 'darkorange'))(100)), cex=0.75, main='sLISA', pch=16)

    return(lisa)
  }))

  ## determine threshold to remove false positives
  mpc <- quantile(lisa, 1-alpha)
  return(mpc)
}
#mpc <- getMinPercentCells(w, m, adjustPv = TRUE, M=10)
mpc <- 0.05

results.filter <- filterSpatialPatterns(mat = m,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = mpc,
                                        verbose = TRUE, details=TRUE)
write.csv(results.filter, file="supp_table_1.csv", quote=FALSE)

results.filter <- filterSpatialPatterns(mat = m,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = mpc,
                                        verbose = TRUE, details=FALSE)
length(results.filter)
length(unlist(dg.sig))
length(intersect(results.filter, unlist(dg.sig)))
length(setdiff(results.filter, unlist(dg.sig)))


## gene that is driven by too few cells
par(mfrow=c(2,2))
#bad.genes <- setdiff(results.filter, unlist(dg.sig))
bad.genes <- setdiff(rownames(I)[I$p.adj < 0.05], results.filter)
length(bad.genes)
g <- bad.genes[1]
gexp <- mat[g,]
moranTest(gexp, w)
hist(gexp)
#sum(lisaTest(gexp, w)$p.value < 0.05)
#plotEmbedding(pos, col=scale(gexp)[,1], cex=2)
#slisa <- MERingue:::getSignedLisa(gexp, w, plot=TRUE, cex=2)
#sum(slisa > -log10(0.05))
par(mfrow=c(5,2), mar=rep(2,4))
invisible(lapply(bad.genes[1:20], function(g) {
  gexp <- scale(mat[g,])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  #invisible(interpolate(pos[colnames(m),], m[g,], main=g, zlim=c(-2,2)))
  plot(pos, col=MERingue:::map2col(gexp), pch=16, cex=0.75, axes=FALSE, frame.plot=TRUE, main=g)
  slisa <- MERingue:::getSignedLisa(gexp, w, plot=FALSE)
  slisa <- scale(slisa)[,1]
  slisa[slisa < -1] <- -1
  slisa[slisa > 1] <- 1
  plot(pos, col=MERingue:::map2col(slisa, pal=colorRampPalette(c('darkgreen', 'white', 'darkorange'))(100)), cex=0.75, main='sLISA', pch=16)
}))

#scc <- spatialCrossCorMatrix(as.matrix(m[bad.genes,]), w)
#m <- t(scale(t(as.matrix(mat))))
#ggroup <- groupSigSpatialPatterns(pos, as.matrix(m[bad.genes,]), scc, power=1, hclustMethod='ward.D', deepSplit = 0)
#gcol <- rainbow(length(levels(ggroup$groups)), v=0.5)[ggroup$groups]
#names(gcol) <- names(ggroup$groups)

## spatial cross correlation
length(results.filter)
scc <- spatialCrossCorMatrix(as.matrix(m[results.filter,]), w)
range(scc)

## double check
plot(I[results.filter,]$observed, diag(scc))
hist(scc, breaks=100)

## group into spatial patterns
#m <- t(scale(t(as.matrix(mat))))
m <- mat
#ggroup <- groupSigSpatialPatterns(pos, as.matrix(m[results.filter,]), scc, power=1, hclustMethod='ward.D2', deepSplit = 2)
ggroup <- groupSigSpatialPatterns(pos, as.matrix(m[results.filter,]), scc, power=1, hclustMethod='ward.D', deepSplit = 4)
gcol <- rainbow(length(levels(ggroup$groups)), v=0.5)[ggroup$groups]
names(gcol) <- names(ggroup$groups)

library(gplots)
bar <- scc
range(bar)
bar[bar < -0.5] <- -0.5
bar[bar > 0.5] <- 0.5
heatmap.2(bar[ggroup$hc$labels,ggroup$hc$labels], trace='none', scale='none',
          Colv=as.dendrogram(ggroup$hc),
          Rowv=as.dendrogram(ggroup$hc),
          labRow=NA, labCol=NA,
          ColSideColors=gcol[ggroup$hc$labels],
          RowSideColors=gcol[ggroup$hc$labels],
          col=colorRampPalette(c('white','black'))(100)
)

#pdf('mOB_genes_all.pdf', width=6, height=10)
par(mfrow=c(5,2), mar=rep(2,4))
invisible(lapply(levels(ggroup$groups), function(x) {
  foo <- names(ggroup$groups)[ggroup$groups==x]
  subfoo <- I[foo,]
  subfoo <- subfoo[order(subfoo$observed, decreasing=TRUE),]
  invisible(lapply(rownames(subfoo)[1:20], function(g) {
    invisible(interpolate(pos[colnames(m),], m[g,], main=g))
  }))
}))
#dev.off()

sp.genes <- unlist(lapply(levels(ggroup$groups), function(x) {
  names(ggroup$groups)[ggroup$groups==x]
}))
library(RColorBrewer)
ccol <- rainbow(length(levels(annot)))[annot]
names(ccol) <- names(annot)
m <- as.matrix(mat[sp.genes,names(sort(annot))])
m <- t(scale(t(m)))
range(m)
m[m < -2.5] <- -2.5
m[m > 2.5] <- 2.5
library(gplots)
ccol <- ccol[colnames(m)]
gcol <- gcol[rownames(m)]
rownames(m)[!(rownames(m) %in% marker.genes)] <- NA
heatmap.2(m, trace="none", Colv=NA, Rowv=NA, labCol=NA,
          ColSideColors=ccol,
          RowSideColors=gcol,
          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)

################ Compare differential expression
diffgexp <- dg.sig
spatgexp <- lapply(levels(ggroup$groups), function(x) {
  names(ggroup$groups)[ggroup$groups==x]
})
names(spatgexp) <- paste0('spatial', levels(ggroup$groups))
spatgexp

library(UpSetR)
upset(UpSetR::fromList(c(diffgexp, spatgexp)), sets=names(c(diffgexp, spatgexp)), keep.order=TRUE, order.by="degree")

## significance of overlap
sigoverlap <- do.call(rbind, lapply(1:length(spatgexp), function(i) {
  so <- unlist(lapply(1:length(diffgexp), function(j) {

    names(spatgexp)[i]
    names(diffgexp)[j]

    #x = # of genes in common between two groups.
    #n = # of genes in group 1.
    #D = # of genes in group 2.
    #N = total genes
    #The representation factor = x / expected # of genes.
    #Expected # of genes = (n * D) / N

    x <- length(intersect(spatgexp[[i]], diffgexp[[j]])) ## shared
    n <- length(spatgexp[[i]])
    D <- length(diffgexp[[j]])
    N <- nrow(counts) ## total
    representation <- x/(n*D/N)

    x
    n
    D
    N
    print(representation)

    phyper(x, D, N-D, n, lower.tail=FALSE)
  }))
  names(so) <- names(diffgexp)
  so
}))
rownames(sigoverlap) <- 1:length(spatgexp)

# visualize
#foo <- p.adjust(sigoverlap)
foo <- sigoverlap
foo <- matrix(foo, nrow(sigoverlap))
rownames(foo) <- paste0('spatial: ', rownames(sigoverlap))
colnames(foo) <- colnames(sigoverlap)
foo[foo < 0.001] <- 0.001
foo <- -log10(foo)
matrix.sort <- function(matrix) {
  row.max <- apply(matrix,1,which.max)
  if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
  matrix[names(sort(row.max)),]
}
foo <- matrix.sort(foo)
heatmap.2(foo, col=colorRampPalette(c('black', 'white', 'red'))(100), scale="none", trace="none", Rowv=NA, Colv=NA, margins = c(25,15))

## plot a few examples
par(mfrow=c(1,1))
plotEmbedding(pos, groups=annot, show.legend=TRUE, cex=3)

par(mfrow=c(5,2), mar=rep(2,4))

g <- intersect(spatgexp[[1]], diffgexp$`1: Granular Cell Layer`)
g <- rownames(I[g[order(I[g,]$p.adj, decreasing=FALSE)],])
g
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))

g <- intersect(spatgexp[[3]], diffgexp$`2: Mitral Cell Layer`)
g <- rownames(I[g[order(I[g,]$p.adj, decreasing=FALSE)],])
g
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))

g <- intersect(spatgexp[[5]], diffgexp$`3: Outer Plexiform Layer`)
g <- rownames(I[g[order(I[g,]$p.adj, decreasing=FALSE)],])
g
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))

g <- intersect(spatgexp[[5]], diffgexp$`4: Glomerular Layer`)
g <- rownames(I[g[order(I[g,]$p.adj, decreasing=FALSE)],])
g
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))

g <- intersect(spatgexp[[4]], diffgexp$`5: Olfactory Nerve Layer`)
g <- rownames(I[g[order(I[g,]$p.adj, decreasing=FALSE)],])
g
invisible(interpolate(pos[colnames(mat),], mat[g[1],], main=g[1]))




################### Induce non-homogeneity

pos_sub <- pos
pos_sub[,1] <- pos_sub[,1]-mean(pos_sub[,1])
vi <- pos_sub[,1] > 0
pos_sub[vi,1] <- 1.1^abs(pos_sub[vi,1])
vi <- pos_sub[,1] < 0
pos_sub[vi,1] <- -1.1^abs(pos_sub[vi,1])
#pos_sub[,1] <- pos_sub[,1]*5
#pos_sub[,1] <- pos_sub[,1]+mean(pos[,1])
plot(pos_sub)

par(mfrow=c(1,1), mar=rep(5,4))
w2 <- voronoiAdjacency(pos_sub, filterDist=2.5, plot=TRUE)
dim(w)
dim(w2)
plotNetwork(pos_sub, w-w2, axes=FALSE, xlab=NA, ylab=NA, main='False Negatives', line.col='green', line.power=2); box() ## difference
plotNetwork(pos_sub, w2-w, axes=FALSE, xlab=NA, ylab=NA, main='False Positives', line.col='green', line.power=2); box() ## difference

I_sub <- getSpatialPatterns(mat, w2)
results.filter_sub <- filterSpatialPatterns(mat = mat,
                                            I = I_sub,
                                            w = w2,
                                            adjustPv = TRUE,
                                            alpha = 0.05,
                                            minPercentCells = mpc,
                                            verbose = TRUE)


length(results.filter_sub)
length(intersect(results.filter_sub, results.filter))
length(setdiff(results.filter_sub, results.filter))
length(setdiff(results.filter, results.filter_sub))

plot(-log10(I$p.value), -log10(I_sub$p.value), pch=16, col=rgb(0,0,0,0.5))
df = cbind(-log10(I$p.value), -log10(I_sub$p.value))
df[is.infinite(df)] <- NA
df <- na.omit(df)
df
ct <- cor.test(df[,1], df[,2])
colnames(df) <- c('x', 'y')
lmfit <- lm(y~x, data.frame(df))
abline(lmfit$coefficients, col='red')
summary(lmfit)
text(5,15,paste0('R2 = ', summary(lmfit)$adj), col='red')
text(5,14,'p-value < 2.2e-16', col='red')

I['Doc2g',]
I_sub['Doc2g',]


library(UpSetR)

length(results.filter_sub)
length(results.filter)
length(intersect(results.filter_sub, results.filter))
length(setdiff(results.filter_sub, results.filter))
length(setdiff(results.filter, results.filter_sub))

foo <- list(results.filter, results.filter_sub)
names(foo) <- c('homogeneous', 'non-homogeneous')
upset(UpSetR::fromList(foo))












################## Compare with overdispersion
library(MUDAN)
varnorm <- normalizeVariance(cd[rownames(mat),], plot=TRUE, details=TRUE)
names(varnorm)
head(varnorm$df)
vg <- varnorm$df$res
names(vg) <- rownames(varnorm$df)
spatg <- I$p.value
names(spatg) <- rownames(I)
plot(-log10(spatg), vg)
#length(spatg)
#length(vg)

################## Spatial heterogeneity within cluster

levels(annot)
sub <- names(annot)[annot == levels(annot)[5]]
w.sub <- getSpatialWeights(pos[sub,], klist=9)
par(mfrow=c(1,1), mar=rep(5,4))
plotNetwork(pos[sub,], w.sub)

gs <- dg.sig[[levels(annot)[5]]]
I.sub <- getSpatialPatterns(mat[gs,sub], w.sub)
results.filter.sub <- filterSpatialPatterns(mat = mat[gs, sub],
                                            I = I.sub,
                                            w = w.sub,
                                            adjustPv = TRUE,
                                            alpha = 0.05,
                                            minPercentCells = 0.05,
                                            verbose = TRUE)
length(results.filter.sub)

g <- results.filter.sub[1]
gexp <- mat[g,sub]
#gexp <- log10(mat[g,]+1)
moranTest(gexp, w.sub)
hist(gexp)
plot(pos[sub,], col=MERingue:::map2col(gexp), pch=16)

## spatial cross correlation
scc.sub <- spatialCrossCorMatrix(as.matrix(mat[results.filter.sub, sub]), w.sub)
ggroup.sub <- groupSigSpatialPatterns(pos[sub,], as.matrix(mat[results.filter.sub, sub]), scc.sub, power = 3, hclustMethod='complete', deepSplit = 1, minClusterSize = 0)

sp.genes.sub <- unlist(lapply(levels(ggroup.sub$groups), function(x) {
  names(ggroup.sub$groups)[ggroup.sub$groups==x]
}))
length(sp.genes.sub)
gcol <- rainbow(length(levels(ggroup.sub$groups)))[ggroup.sub$groups]
names(gcol) <- names(ggroup.sub$groups)
m <- as.matrix(mat[sp.genes.sub,sub])
m <- t(scale(t(m)))
range(m)
m[m < -2.5] <- -2.5
m[m > 2.5] <- 2.5
library(gplots)
heatmap.2(m, trace="none", Colv=NA, Rowv=NA, labRow=NA, labCol=NA,
          RowSideColors=gcol[rownames(m)],
          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)
heatmap.2(m, trace="none", Rowv=NA, labRow=NA, labCol=NA,
          RowSideColors=gcol[rownames(m)],
          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)

library(gplots)
heatmap.2(scc.sub[sp.genes.sub, sp.genes.sub], trace="none", Colv=NA, Rowv=NA, labRow=NA, labCol=NA,
          ColSideColors=gcol[rownames(m)],
          RowSideColors=gcol[rownames(m)],
          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
          scale='none'
)

length(sp.genes.sub)
length(sp.genes)
length(setdiff(sp.genes.sub, sp.genes))
ggroup$groups[sp.genes.sub]

I[setdiff(sp.genes.sub, sp.genes),]
g <- "Vps13c"
g <- "Gas2l3"
I.sub[g,]
I[g,]

gexp <- mat[g,]
plot(pos, col=MERingue:::map2col(gexp), pch=16)
gexp <- mat[g,sub]
plot(pos[sub,], col=MERingue:::map2col(gexp), pch=16)

spatgexp <- lapply(levels(ggroup$groups), function(x) {
  names(ggroup$groups)[ggroup$groups==x]
})
names(spatgexp) <- paste0('spatial', levels(ggroup$groups))
spatgexp.sub <- lapply(levels(ggroup.sub$groups), function(x) {
  names(ggroup.sub$groups)[ggroup.sub$groups==x]
})
names(spatgexp.sub) <- paste0('spatial.sub', levels(ggroup.sub$groups))
upset(UpSetR::fromList(c(spatgexp, spatgexp.sub)), sets=names(c(spatgexp, spatgexp.sub)), keep.order=TRUE, order.by="degree")
