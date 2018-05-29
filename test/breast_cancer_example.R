# Parse data
parseData <- function(jsfile) {
  library(jsonlite)
  result <- jsonlite::fromJSON(jsfile, flatten=TRUE)

  pos <- data.frame(result$barcode, x=result$x, y=result$y)
  vi <- !duplicated(pos)
  table(vi)
  pos <- pos[vi,]
  rownames(pos) <- pos[,1]
  pos <- pos[,-1]

  dat <- data.frame(cell=result$barcode, gene=result$gene, value=result$hits)
  library(reshape2)
  cd <- dcast(data = dat, formula = cell~gene, fun.aggregate = sum, value.var = "value")
  rownames(cd) <- cd[,1]
  cd <- cd[,-1]
  cd <- t(cd)



  return(list(cd=cd, pos=pos))
}
BCL1 <- parseData('Data/BC_layer1.json_.gz')
BCL2 <- parseData('Data/BC_layer2.json_.gz')
BCL3 <- parseData('Data/BC_layer3.json_.gz')
BCL4 <- parseData('Data/BC_layer4.json_.gz')

# Get highly spatially variable genes (union)
getI <- function(BCL1) {
  w <- getSpatialWeights(BCL1$pos, klist=6)
  plotNetwork(pos=BCL1$pos, adj=w)
  I <- getSpatialPatterns(BCL1$cd, w)
  results <- I
  vi <- results$p.adj < 0.05
  vi[is.na(vi)] <- FALSE
  results.sig <- rownames(results)[vi]
  return(results.sig)
}
sig1 <- getI(BCL1)
sig2 <- getI(BCL2)
sig3 <- getI(BCL3)
sig4 <- getI(BCL4)
sig.genes <- unique(c(sig1, sig2, sig3, sig4))
sig.genes

