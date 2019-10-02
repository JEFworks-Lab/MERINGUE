#' Differential expression analysis (adapted from PAGODA2)
#'
#' @param cd A read count matrix. The rows correspond to genes, columns correspond to individual cells
#' @param cols Column/cell group annotations. Will perform one vs. all differential expression analysis.
#' @param verbose Verbosity
#'
#' @export
#'
getDifferentialGenes <- function(cd, cols, verbose=TRUE) {
  cm <- t(cd)

  ## match matrix rownames (cells) and group annotations
  if(!all(rownames(cm) %in% names(cols))) { warning("Cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")}
  ## determine a subset of cells that's in the cols and cols[cell]!=NA
  valid.cells <- rownames(cm) %in% names(cols)[!is.na(cols)];
  if(!all(valid.cells)) {
    ## take a subset of the count matrix
    cm <- cm[valid.cells,]
  }
  ## reorder cols
  cols <- as.factor(cols[match(rownames(cm),names(cols))]);
  cols <- as.factor(cols);

  if(verbose) {
    print(paste0("Running differential expression with ",length(levels(cols))," clusters ... "))
  }

  ## modified from pagoda2
  ## run wilcoxon test comparing each group with the rest
  ## calculate rank per column (per-gene) average rank matrix
  xr <- apply(cm, 2, function(foo) {
    #foo[foo==0] <- NA
    bar <- rank(foo)
    #bar[is.na(foo)] <- 0
    bar[foo==0] <- 0
    bar
  }); rownames(xr) <- rownames(cm)
  # xr <- sparse_matrix_column_ranks(cm);

  ## calculate rank sums per group
  grs <- do.call(rbind, lapply(levels(cols), function(g) Matrix::colSums(xr[cols==g,])))
  rownames(grs) <- levels(cols); colnames(grs) <- colnames(xr)
  # grs <- colSumByFac(xr,as.integer(cols))[-1,,drop=F]

  ## calculate number of non-zero entries per group
  gnzz <- do.call(rbind, lapply(levels(cols), function(g) Matrix::colSums(xr[cols==g,]>0)))
  rownames(gnzz) <- levels(cols); colnames(gnzz) <- colnames(xr)
  # xr@x <- numeric(length(xr@x))+1
  # gnzz <- colSumByFac(xr,as.integer(cols))[-1,,drop=F]

  # group.size <- as.numeric(tapply(cols,cols,length));
  # group.size <- as.numeric(tapply(cols,cols,length))[1:nrow(gnzz)]; group.size[is.na(group.size)]<-0; # trailing empty levels are cut off by colSumByFac
  group.size <- as.numeric(table(cols))

  # add contribution of zero entries to the grs
  gnz <- (group.size-gnzz)

  ## rank of a 0 entry for each gene
  # zero.ranks <- (nrow(xr)-diff(xr@p)+1)/2 # number of total zero entries per gene
  zero.ranks <- apply(cm, 2, function(foo) {
    bar <- rank(foo)
    bar[foo==0][1]
  })
  # if nothing ranked 0, will be NA so fix
  zero.ranks[is.na(zero.ranks)] <- 0

  ustat <- t((t(gnz)*zero.ranks)) + grs - group.size*(group.size+1)/2

  # standardize
  n1n2 <- group.size*(nrow(cm)-group.size);
  # usigma <- sqrt(n1n2*(nrow(cm)+1)/12) # without tie correction
  ## correcting for 0 ties, of which there are plenty
  # usigma <- sqrt(n1n2*(nrow(cm)+1)/12)
  usigma <- sqrt((nrow(cm) +1 - (gnz^3 - gnz)/(nrow(cm)*(nrow(cm)-1)))*n1n2/12)
  # standardized U value- z score
  x <- t((ustat - n1n2/2)/usigma);

  # correct for multiple hypothesis
  y <- matrix(bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE,
                              log.p = TRUE), log = TRUE), ncol = ncol(x)) * sign(x)
  #y <- matrix(pnorm(as.numeric(abs(x)), lower.tail = FALSE,
  #                            log.p = TRUE), ncol = ncol(x)) * sign(x)
  y <- exp(y) # log p-values are difficult to interpret
  y[y>1] <- 1
  x <- matrix(qnorm(bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE,
                                    log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE),
              ncol = ncol(x)) * sign(x)
  rownames(y) <- rownames(x) <- colnames(cm)
  colnames(y) <- colnames(x) <- levels(cols)[1:ncol(x)]

  # add fold change information
  log.gene.av <- log2(Matrix::colMeans(cm));
  group.gene.av <- do.call(rbind, lapply(levels(cols), function(g) Matrix::colSums(cm[cols==g,]>0))) / (group.size+1);
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av;
  # fraction of cells expressing
  f.expressing <- t(gnzz / group.size);
  max.group <- max.col(log2.fold.change)

  if(verbose) {
    print("Summarizing results ... ")
  }

  ## summarize
  ds <- lapply(1:ncol(x),function(i) {
    r <- data.frame(p.adj=y[,i],Z=x[,i],M=log2.fold.change[,i],highest=max.group==i,fe=f.expressing[,i])
    #r <- data.frame(p.value=y[,i],Z=x[,i],M=log2.fold.change[,i],highest=max.group==i,fe=f.expressing[,i])
    rownames(r) <- rownames(x)
    r
  })
  names(ds)<-colnames(x)

  return(ds)
}
## BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x <- x[nai]
  id <- order(x, decreasing = FALSE)
  if(log) {
    q <- x[id] + log(length(x)/seq_along(x))
  } else {
    q <- x[id]*length(x)/seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai] <- a
  ox
}


#' Winsorize expression values to prevent outliers
#'
#' @param x Values
#' @param qt Values below this quantile and above 1-this quantile will be set to the quantile value
#'
#' @export
#'
winsorize <- function (x, qt=.05) {
  if(length(qt) != 1 || qt < 0 ||
     qt > 0.5) {
    stop("bad value for quantile threashold")
  }
  lim <- quantile(x, probs=c(qt, 1-qt))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}
