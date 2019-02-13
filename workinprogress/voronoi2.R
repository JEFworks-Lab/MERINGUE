# get neighbor's neighbor
set.seed(0)
library(MERingue)
data(mOB)

# Clean and normalize data
pos <- mOB$pos
cd <- mOB$counts

counts <- cleanCounts(cd, min.reads=10, min.lib.size=10)
pos <- pos[colnames(counts),]
mat <- normalizeCounts(counts, log=FALSE)

voronoiAdjacency2 <- function(pos, distMethod = "euc", filterDist = NA, plot=FALSE){

  makeixy <- function(data, formula, scale){
    m = model.frame(formula, data=data)
    if(ncol(m)!=3){
      stop("incorrect adjacency formula: id~x+y needed")
    }
    names(m)=c("id","x","y")
    m[,2]=m[,2]/scale
    m[,3]=m[,3]/scale
    m
  }

  data <- data.frame('i'=1:nrow(pos), pos)
  data <- makeixy(data, formula=i~x+y, scale=1)

  P <- nrow(data);  # number of rows

  dd1 = deldir::deldir(data$x,
                       data$y,
                       suppressMsge=TRUE, plotit = FALSE, sort=FALSE);  # find adjacencies
  ## create distance matrix
  D1 <- matrix(0,P,P);
  #D1[as.matrix(dd1$delsgs[,c("ind1","ind2")])] = sqrt((dd1$delsgs[,c("x1")]-dd1$delsgs[,c("x2")])^2+(dd1$delsgs[,c("y1")]-dd1$delsgs[,c("y2")])^2);
  #D1[as.matrix(dd1$delsgs[,c("ind2","ind1")])] = sqrt((dd1$delsgs[,c("x1")]-dd1$delsgs[,c("x2")])^2+(dd1$delsgs[,c("y1")]-dd1$delsgs[,c("y2")])^2);
  a = dd1$delsgs[, c("x1", "y1")]
  b = dd1$delsgs[, c("x2", "y2")]
  colnames(a) <- colnames(b) <- c('x', 'y')
  ds <- sapply(seq_len(nrow(a)), function(i) dist(rbind(a[i,], b[i,]), method=distMethod))
  D1[as.matrix(dd1$delsgs[,c("ind1","ind2")])] <- ds
  D1[as.matrix(dd1$delsgs[,c("ind2","ind1")])] <- ds
  D <- D1

  if(!is.na(filterDist)) {
    D[D>filterDist] = 0
  }

  # binarize
  D[D>0] <- 1

  if(plot) {
    plotNetwork(pos, adj=D)
  }

  rownames(D) <- rownames(pos)
  colnames(D) <- rownames(pos)

  return(D);
}
par(mfrow=c(1,1))
w <- voronoiAdjacency2(pos, distMethod='euc', filterDist = 2, plot=TRUE)
w <- voronoiAdjacency2(pos, distMethod='man', filterDist = 3, plot=TRUE)
w <- voronoiAdjacency2(pos, distMethod='max', filterDist = 3, plot=TRUE)
w <- voronoiAdjacency2(pos, distMethod='min', filterDist = NA, plot=TRUE)
