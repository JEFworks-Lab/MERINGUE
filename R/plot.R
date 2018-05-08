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
plotEmbedding <- function(emb, groups=NULL, colors=NULL, cex=0.6, alpha=0.4, gradientPalette=NULL, zlim=NULL, s=1, v=0.8, min.group.size=1, show.legend=FALSE, mark.clusters=FALSE, mark.cluster.cex=2, shuffle.colors=F, legend.x='topright', gradient.range.quantile=0.95, verbose=TRUE, unclassified.cell.color='gray70', group.level.colors=NULL, ...) {

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

  plot(emb,col=adjustcolor(cols,alpha.f=alpha),cex=cex,pch=19,axes=F, ...); box();
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

