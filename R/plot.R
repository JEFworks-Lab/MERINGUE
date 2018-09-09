#' @import stats grDevices graphics
NULL

#' Plot 2D embedding
#'
#' @param emb dataframe with x and y coordinates
#' @param groups factor annotations for rows on emb for visualizing cluster annotations
#' @param colors color or numeric values for rows on emb for visualizing gene expression
#' @param cex point size
#' @param alpha point opacity
#' @param gradientPalette palette for colors if numeric values provided
#' @param zlim range for colors
#' @param s saturation of rainbow for group colors
#' @param v value of rainbow for group colors
#' @param min.group.size minimum size of group in order for group to be colored
#' @param show.legend whether to show legend
#' @param mark.clusters whether to mark clusters with name of cluster
#' @param mark.cluster.cex cluster marker point size
#' @param shuffle.colors whether to shuffle group colors
#' @param legend.x legend position ie. 'topright', 'topleft', 'bottomleft', 'bottomright'
#' @param gradient.range.quantile quantile for mapping colors to gradient palette
#' @param verbose verbosity
#' @param unclassified.cell.color cells not included in groups will be labeled in this color
#' @param group.level.colors set group level colors. Default uses rainbow.
#' @param ... Additional parameters to pass to BASE::plot
#'
#' @export
plotEmbedding <- function(emb, groups=NULL, colors=NULL, cex=0.6, alpha=0.4, gradientPalette=NULL, zlim=NULL, s=1, v=0.8, min.group.size=1, show.legend=FALSE, mark.clusters=FALSE, mark.cluster.cex=2, shuffle.colors=F, legend.x='topright', gradient.range.quantile=0.95, verbose=TRUE, unclassified.cell.color='gray70', group.level.colors=NULL, xlab=NA, ylab=NA, ...) {

  if(!is.null(colors)) {
    ## use clusters information
    if(!all(rownames(emb) %in% names(colors))) { warning("provided cluster vector doesn't list colors for all of the cells; unmatched cells will be shown in gray. ")}
    if(all(areColors(colors))) {
      if(verbose) cat("using supplied colors as is\n")
      cols <- colors[match(rownames(emb),names(colors))]; cols[is.na(cols)] <- unclassified.cell.color;
      names(cols) <- rownames(emb)
    } else {
      if(is.numeric(colors)) { # treat as a gradient
        if(verbose) cat("treating colors as a gradient")
        if(is.null(gradientPalette)) { # set up default gradients
          if(all(sign(colors)>=0)) {
            gradientPalette <- colorRampPalette(c('gray80','red'), space = "Lab")(1024)
          } else {
            gradientPalette <- colorRampPalette(c("blue", "grey70", "red"), space = "Lab")(1024)
          }
        }
        if(is.null(zlim)) { # set up value limits
          if(all(sign(colors)>=0)) {
            zlim <- as.numeric(quantile(colors,p=c(1-gradient.range.quantile,gradient.range.quantile)))
            if(diff(zlim)==0) {
              zlim <- as.numeric(range(colors))
            }
          } else {
            zlim <- c(-1,1)*as.numeric(quantile(abs(colors),p=gradient.range.quantile))
            if(diff(zlim)==0) {
              zlim <- c(-1,1)*as.numeric(max(abs(colors)))
            }
          }
        }
        # restrict the values
        colors[colors<zlim[1]] <- zlim[1]; colors[colors>zlim[2]] <- zlim[2];

        if(verbose) cat(' with zlim:',zlim,'\n')
        colors <- (colors-zlim[1])/(zlim[2]-zlim[1])
        cols <- gradientPalette[colors[match(rownames(emb),names(colors))]*(length(gradientPalette)-1)+1]
        names(cols) <- rownames(emb)
      } else {
        stop("colors argument must be a cell-named vector of either character colors or numeric values to be mapped to a gradient")
      }
    }
  } else {
    if(!is.null(groups)) {
      if(min.group.size>1) { groups[groups %in% levels(groups)[unlist(tapply(groups,groups,length))<min.group.size]] <- NA; groups <- droplevels(groups); }
      groups <- as.factor(groups)[rownames(emb)]
      if(verbose) cat("using provided groups as a factor\n")
      factor.mapping=TRUE;
      ## set up a rainbow color on the factor
      factor.colors <- fac2col(groups,s=s,v=v,shuffle=shuffle.colors,min.group.size=min.group.size,unclassified.cell.color=unclassified.cell.color,level.colors=group.level.colors,return.details=T)
      cols <- factor.colors$colors;
      names(cols) <- rownames(emb)
    } else {
      cols <- rep(unclassified.cell.color, nrow(emb))
      names(cols) <- rownames(emb)
    }
  }

  plot(emb,col=adjustcolor(cols,alpha.f=alpha),cex=cex,pch=19,axes=F,xlab=xlab,ylab=ylab, ...); box();
  if(mark.clusters) {
    if(!is.null(groups)) {
      cent.pos <- do.call(rbind,tapply(1:nrow(emb),groups,function(ii) apply(emb[ii,,drop=F],2,median)))
      cent.pos <- na.omit(cent.pos);
      text(cent.pos[,1],cent.pos[,2],labels=rownames(cent.pos),cex=mark.cluster.cex)
    }
  }
  if(show.legend) {
    if(factor.mapping) {
      legend(x=legend.x,pch=rep(19,length(levels(groups))),bty='n',col=factor.colors$palette,legend=names(factor.colors$palette))
    }
  }
}

## a utility function to translate factor into colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.details=F,unclassified.cell.color='gray50',level.colors=NULL) {
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
    x <- droplevels(x)
  }
  if(is.null(level.colors)) {
    col <- rainbow(length(levels(x)),s=s,v=v);
  } else {
    col <- level.colors[1:length(levels(x))];
  }
  names(col) <- levels(x);

  if(shuffle) col <- sample(col);

  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  if(return.details) {
    return(list(colors=y,palette=col))
  } else {
    return(y);
  }
}

## quick utility to check if given character vector is colors
## thanks, Josh O'Brien: http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
areColors <- function(x) {
  is.character(x) &
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
  })
}








#' Plot neighbor network
#' https://stackoverflow.com/questions/43879347/plotting-a-adjacency-matrix-using-pure-r
#'
#' @export
plotNetwork <- function(pos, adj, col='black', line.col='red', line.power=1, ...) {
  plot(pos, pch=16, col=col, ...)
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


#' Helper function to map values to colors
#' Source: https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
map2col <- function(x, pal=colorRampPalette(c('blue', 'grey', 'red'))(100), na.col='lightgrey', limits=NULL){
  original <- x
  x <- na.omit(x)
  if(is.null(limits)) limits=range(x)
  y <- pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  names(y) <- names(x)

  colors <- rep(na.col, length(original))
  names(colors) <- names(original)
  colors[names(y)] <- y

  return(colors)
}



#' Gridded bivariate interpolation
#' @export
interpolate <- function(pos, gexp, trim=0, zlim=c(-1,1), fill=TRUE, binSize=100, cex=1, col=colorRampPalette(c('blue', 'white', 'red'))(100), plot=TRUE, ...) {
    z <- scale(winsorize(gexp, trim))[,1]
    names(z) <- names(gexp)
    z[z < zlim[1]] <- zlim[1]
    z[z > zlim[2]] <- zlim[2]
    x <- pos[,1]
    y <- pos[,2]
    names(x) <- names(y) <- rownames(pos)

    if(fill) {
        zb <- rep(0, nrow(pos))
        names(zb) <- rownames(pos)
        zb[names(gexp)] <- z
    } else {
        x <- x[names(gexp)]
        y <- y[names(gexp)]
        zb <- z
    }

    int <- akima::interp(x, y, zb, nx=binSize, ny=binSize, linear=TRUE)

    if(plot) {
        plot(pos[names(gexp),], col=map2col(z), pch=16, axes=FALSE, xlab=NA, ylab=NA, ...); box()
        image(int, col=col, axes=FALSE, frame.plot=TRUE, ...)
    }

    return(int)
}


#' 2D kernel density estimation
plotInterpolatedDensity <- function(pos) {
  x <- pos[,1]
  y <- pos[,2]
  dens <- MASS::kde2d(x, y)
  persp(dens, phi = 30, theta = 20, d = 5)
}



#' Get density of points in 2 dimensions.
#' @author Kamil Slowikowski http://slowkow.com/notes/ggplot2-color-by-density/
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @return The density within each square.
#'
getDensity <- function(pos, n = 100) {
  x <- pos[,1]
  y <- pos[,2]
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  dens <- dens$z[ii]
  names(dens) <- rownames(pos)
  return(dens)
}

plotDensity <- function(pos, n=100) {
  dens <- getDensity(pos, n)
  plot(pos, col=map2col(dens), pch=16)
}

plotDensitySubset <- function(pos, subset, n=100) {
  dens <- getDensity(pos[subset,], n)
  densall <- rep(NA, nrow(pos))
  names(densall) <- rownames(pos)
  densall[subset] <- dens
  plot(pos, col=map2col(densall), pch=16)
}




#' plot boxplot of expression for nearest neighbor
#' @export
plotNeighborBoxplot <- function(gexpA, gexpB, groupA, groupB, weight) {
    par(mfrow=c(1,2))

    ## plot correlation between groupA cells and neighbors
    nbs <- lapply(groupA, function(x) names(which(weight[x,]==1)))
    names(nbs) <- groupA
    ## gene A expression in group A
    foo <- gexpA[groupA]

    ## average gene B expression for neighbors from group B
    bar <- unlist(lapply(nbs, function(y) mean(gexpB[y])))

    lmfit <- lm(bar~foo)

    plot(foo, bar, xlab='gexpA in groupA', ylab='gexpB in nearest neighbor in groupB')
    abline(lmfit)
    pv <- coef(summary(lmfit))['foo','Pr(>|t|)']
    text(0,0, paste0('p-value: ', signif(pv, 3)), adj=0, col='red', font=2)

    ## boxplot of distribution
    bard <- do.call(rbind, lapply(names(nbs), function(y) {
                               ngexp <- gexpB[nbs[[y]]]
                               cbind(cell=rep(y, length(ngexp)), ngexp=ngexp)
                           }))
    bard <- cbind(bard, gexp=foo[bard[,'cell']])
    bard <- data.frame(bard)
    class(bard$ngexp) <- 'numeric'
    class(bard$gexp) <- 'numeric'
    boxplot(ngexp~gexp, data=bard)

    return(lmfit)
}




#' Rotate position by angle theta in radians
rotatePos <- function(pos, theta) {
  rotMat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2)
  pos2 <- t(rotMat %*% t(pos))
  colnames(pos2) <- colnames(pos)
  return(pos2)
}
