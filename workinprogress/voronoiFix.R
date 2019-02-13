#' Neighbor weight matrix by voronoi adjacency
#'
#' @description Neighbor weight matrix by voronoi adjacency. Modified from Barry Rowlingson and Alison Hale's voronoi_adjacency function in the caramellar package
#'
#' @param pos Position
#' @param filterDist Euclidean distance beyond which two cells cannot be considered neighbors
#' @param plot Whether to plot
#'
#' @examples {
#' data(mOB)
#' pos <- mOB$pos
#' w <- voronoiAdjacency(pos, plot=TRUE)
#' }
#'
#' @export
#'
voronoiAdjacency <- function(pos, filterDist = NA, plot=FALSE){

  if(ncol(pos)>2) {
    stop('2D tesselation only')
  }
  if(is.null(colnames(pos))) {
    colnames(pos) <- c('x', 'y')
  }
  if(is.null(rownames(pos))) {
    rownames(pos) <- seq_len(nrow(pos))
  }

  ## deldir doesn't handle edge cases well
  ## create dummy edges
  pos0 <- pos
  t <- max((max(pos[,1])-min(pos[,1]))/10, (max(pos[,2])-min(pos[,2]))/10)
  x1 <- min(pos[,1])-t
  x2 <- max(pos[,1])+t
  y1 <- min(pos[,2])-t
  y2 <- max(pos[,2])+t
  dummy <- rbind(
    cbind(seq(x1, x2, by=(x2-x1)/10), y1),
    cbind(seq(x1, x2, by=(x2-x1)/10), y2),
    cbind(x1, seq(y1, y2, by=(y2-y1)/10)),
    cbind(x2, seq(y1, y2, by=(y2-y1)/10))
  )
  colnames(dummy) <- c('x', 'y')
  rownames(dummy) <- paste0('dummy', 1:nrow(dummy))
  pos <- rbind(pos, dummy)

  ## find adjacencies
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
  data <- data.frame('i'=seq_len(nrow(pos)), pos)
  data <- makeixy(data, formula=i~x+y, scale=1)
  dd1 = deldir::deldir(data$x,
                       data$y,
                       suppressMsge=TRUE, sort=FALSE,
                       plotit = plot, col='grey')
  if(plot) { box() }

  ## create distance matrix
  P <- nrow(data) # number of rows
  D1 <- matrix(0,P,P)
  rownames(D1) <- rownames(pos)
  colnames(D1) <- rownames(pos)
  D1[as.matrix(dd1$delsgs[,c("ind1","ind2")])] = sqrt((dd1$delsgs[,c("x1")]-dd1$delsgs[,c("x2")])^2+(dd1$delsgs[,c("y1")]-dd1$delsgs[,c("y2")])^2)
  D1[as.matrix(dd1$delsgs[,c("ind2","ind1")])] = sqrt((dd1$delsgs[,c("x1")]-dd1$delsgs[,c("x2")])^2+(dd1$delsgs[,c("y1")]-dd1$delsgs[,c("y2")])^2)
  D <- D1

  ## filter by distance
  if(!is.na(filterDist)) {
    D[D>filterDist] = 0
  }
  ## filter to original set of points
  D <- D[rownames(pos0), rownames(pos0)]

  ## binarize
  D[D>0] <- 1

  ## plot
  if(plot) {
    idx <- which(D>0, arr.ind = T)
    for(i in seq_len(nrow(idx))) {
      lines(
        c(pos0[idx[i,1],1], pos0[idx[i,2],1]),
        c(pos0[idx[i,1],2], pos0[idx[i,2],2]),
        col='red'
      )
    }
    points(pos0, pch=16)
  }


  return(D);
}
