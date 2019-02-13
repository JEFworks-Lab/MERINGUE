## sup 4

############# LISA

par(mfrow=c(1,1))

set.seed(0)
N <- 10^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
plotEmbedding(pos)

pos[,1] <- expPos(pos[,1])
pos[,2] <- expPos(pos[,2])
plot(pos)

set.seed(0)
posj <- jitter(pos, amount=0.015)
plot(posj)

#gexp <- scale(posj[,1])[,1]

#gexp <- rnorm(N)
gexp <- rep(-1, N)
vi <- posj[,1] > 0.8 & posj[,2] > 0.8
gexp[vi] <- 1

# noise
vi <- posj[,1] < 0.4 & posj[,2] < 0.4
gexp[vi] <- 1

#gexp <- scale(gexp)[,1]
range(gexp)
#gexp <- jitter(gexp, amount=1)
names(gexp) <- rownames(posj)

w <- voronoiAdjacency(posj)

li <- lisaTest(gexp, w)
li
lipv <- -log10(li$p.value)*sign(gexp)
lipv[is.infinite(lipv)] <- NA
lipv[is.na(lipv)] <- max(lipv, na.rm=TRUE)
range(lipv)
lipv[lipv < -max(lipv)] <- -max(lipv)
lipv[lipv > -min(lipv)] <- -min(lipv)
#lipv[is.infinite(lipv)] <- 20
#lipv <- scale(lipv)[,1]
names(lipv) <- names(gexp)
range(lipv)

par(mfrow=c(1,2))
plotEmbedding(posj, col=gexp, cex=2, main = 'Gene Expression')
plotEmbedding(posj, col=lipv, cex=2, main='signed LISA')



############## Correlogram

voronoiAdjacency <- function(pos, filterDist = NA, plot=FALSE, epsilon=5){

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
  t <- max((max(pos[,1])-min(pos[,1]))/epsilon, (max(pos[,2])-min(pos[,2]))/epsilon)
  x1 <- min(pos[,1])-t
  x2 <- max(pos[,1])+t
  y1 <- min(pos[,2])-t
  y2 <- max(pos[,2])+t
  dummy <- rbind(
    cbind(seq(x1, x2, by=(x2-x1)/epsilon), y1),
    cbind(seq(x1, x2, by=(x2-x1)/epsilon), y2),
    cbind(x1, seq(y1, y2, by=(y2-y1)/epsilon)),
    cbind(x2, seq(y1, y2, by=(y2-y1)/epsilon))
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


correlogram <- function(gexp, w, posj, k=sqrt(nrow(posj)), t=0.05) {

  correlogramHelper <- function(gexp, w) {
    #plotEmbedding(posj, col=gexp, cex=1)
    #plotNetwork(posj, w)

    diag(w) <- 0
    w2 <- matrix(0, nrow(w), ncol(w))
    rownames(w2) <- colnames(w2) <- rownames(w)

    for(x in seq_along(rownames(w))) {
      delta <- w[x,]==1 # neighbors
      if(sum(delta) > 1) {
        echo <- w[delta,]==1 # neighbor's neighbors
        if(nrow(echo)>1) {
          foxtrot <- colSums(echo)>0
        } else {
          foxtrot <- echo>0
        }
        w2[x, foxtrot] <- 1
      }
    }

    # neighbors of neighbors driving pattern
    #li <- lisaTest(gexp, w)$p.value
    #vi <- li > t
    #names(vi) <- names(gexp)
    #plotEmbedding(posj, groups=vi)

    #weight <- w
    #plotNetwork(posj, weight)
    #weight[!vi,!vi] <- 0
    #plotNetwork(posj, weight)

    diag(w2) <- 0
    #w2[w==1] <- 0
    obs <- MERingue:::moranTest(gexp, w2)$obs # observed
    return(list(obs=obs, w2=w2))
  }

  # init
  w <- voronoiAdjacency(posj)
  x <- MERingue:::moranTest_C(gexp, w)[1]

  # neighbors of neighbors
  #par(mfrow=c(k,3))
  for(i in seq_len(k)) {

    #plotEmbedding(posj, col=gexp, cex=2, main = 'Gene Expression')
    #plotNetwork(posj, w)

    # who are their neighbors
    results <- correlogramHelper(gexp, w)
    x <- c(x, results$obs)
    x
    w <- results$w2

    # neighbors of neighbors driving pattern
    #li <- lisaTest(gexp, w)$p.value
    #vi <- li > t
    #names(vi) <- names(gexp)
    #plotEmbedding(posj, groups=vi, cex=1)

  }

  return(x)
}


correlogram2 <- function(gexp, w, posj, k=sqrt(nrow(posj))/2, t=0.05) {
  #plot(posj)
  w <- voronoiAdjacency(posj)
  #diag(w)
  #plotNetwork(posj, w)

  # driving pattern
  li <- lisaTest(gexp, w)$p.value
  vi <- li <= t
  names(vi) <- names(gexp)
  drivers <- names(vi)[which(vi)]
  #plotEmbedding(posj, groups=vi)

  correlogram2Helper <- function(w) {
    w2 <- do.call(rbind, lapply(rownames(w), function(i) {
      #print(i)
      foo <- posj
      vi <- w[i, ]!=1
      foo <- foo[vi,]
      #plot(foo)
      w2 <- voronoiAdjacency(foo)
      #plotNetwork(foo, w2)
      x <- w2[i,]
      y <- x[rownames(w)]
      names(y) <- rownames(w)
      y[is.na(y)] <- 0

      return(y)
    }))
    rownames(w2) <- rownames(w)
    #plotNetwork(posj, foo)

    #obs
    return(w2)
  }

  trimNetwork <- function(w, drivers) {
    for(i in seq_len(nrow(w))) {
      w[i, setdiff(names(which(w[i,]==1)), drivers)] <- 0
    }
    return(w)
    #plotNetwork(posj, w)
  }

  # init
  w <- voronoiAdjacency(posj)
  #plotNetwork(posj, w)
  x <- MERingue:::moranTest_C(gexp, w)[1]
  # neighbors of neighbors
  for(i in seq_len(k)) {
    print(i)
    w2 <- correlogram2Helper(w)
    #plotNetwork(posj, w2)
    w0 <- trimNetwork(w2, drivers)
    plotNetwork(posj, w0, main=i)
    drivers <- c(drivers, unique(rownames(which(w0==1, arr.ind=TRUE))))
    obs <- MERingue:::moranTest(gexp, w0)$obs # observed
    x <- c(x, obs)
    w <- w2
  }
  return(x)
}



# check
par(mfrow=c(1,1))
set.seed(0)
N <- 10^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
plotEmbedding(pos)

#pos[,1] <- expPos(pos[,1])
#pos[,2] <- expPos(pos[,2])
#plot(pos)

set.seed(0)
posj <- jitter(pos, amount=0.015)
plot(posj)

gexp <- rep(-1, N)
vi <- posj[,1] > max(pos[,1])*0.7 & posj[,2] > max(pos[,2])*0.7
gexp[vi] <- 1
names(gexp) <- rownames(posj)
x <- correlogram2(gexp, w, posj)

gexp2 <- scale(posj[,1]*posj[,2], center=TRUE, scale=FALSE)[,1]
names(gexp2) <- rownames(posj)
x2 <- correlogram2(gexp2, w, posj)

par(mfrow=c(2,2))
plotEmbedding(posj, color=gexp)
plot(1:length(x), x, type="l")
plotEmbedding(posj, color=gexp2)
plot(1:length(x2), x2, type="l")


