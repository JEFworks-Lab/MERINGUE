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


#' Get overdispersed genes whose expression variance is higher than transcriptome-wide expectations
#' (Modified from SCDE/PAGODA2 code)
#'
#' Normalizes gene expression magnitudes to with respect to its ratio to the
#' transcriptome-wide expectation as determined by local regression on expression magnitude
#'
#' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
#' @param gam.k Generalized additive model parameter; the dimension of the basis used to represent the smooth term (default: 5)
#' @param alpha Significance threshold (default: 0.05)
#' @param plot Whether to plot the results (default: FALSE)
#' @param use.unadjusted.pvals If true, will apply BH correction (default: FALSE)
#' @param do.par Whether to adjust par for plotting if plotting (default: TRUE)
#' @param max.adjusted.variance Ceiling on maximum variance after normalization to prevent infinites (default: 1e3)
#' @param min.adjusted.variance Floor on minimum variance after normalization (default: 1e-3)
#' @param verbose Verbosity (default: TRUE)
#' @param details If true, will return data frame of normalization parameters. Else will return normalized matrix.(default: FALSE)
#'
#' @return If details is true, will return data frame of normalization parameters. Else will return overdispersed genes.
#'
#' @examples
#' data(mOB)
#' ods <- getOverdispersedGenes(mOB$counts)
#'
#' @importFrom mgcv s
#'
#' @export
#'
getOverdispersedGenes <- function(counts, gam.k=5, alpha=0.05, plot=FALSE, use.unadjusted.pvals=FALSE, do.par=TRUE, max.adjusted.variance=1e3, min.adjusted.variance=1e-3, verbose=TRUE, details=FALSE) {

  if (!class(counts) %in% c("dgCMatrix", "dgTMatrix")) {
    if(verbose) {
      message('Converting to sparse matrix ...')
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  mat <- Matrix::t(counts) ## make rows as cells, cols as genes

  if(verbose) {
    print("Calculating variance fit ...")
  }
  dfm <- log(Matrix::colMeans(mat))
  dfv <- log(apply(mat, 2, var))
  names(dfm) <- names(dfv) <- colnames(mat)
  df <- data.frame(m=dfm, v=dfv)

  vi <- which(is.finite(dfv))

  if(length(vi)<gam.k*1.5) { gam.k=1 } ## too few genes

  if(gam.k<2) {
    if(verbose) {
      print("Using lm ...")
    }
    m <- lm(v ~ m, data = df[vi,])
  } else {
    if(verbose) {
      print(paste0("Using gam with k=", gam.k, "..."))
    }
    fm <- as.formula(sprintf("v ~ s(m, k = %s)", gam.k))
    m <- mgcv::gam(fm, data = df[vi,])
  }
  df$res <- -Inf;  df$res[vi] <- resid(m,type='response')
  n.cells <- ncol(mat)
  n.obs <- nrow(mat)
  df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=T))
  df$lpa <- bh.adjust(df$lp,log=TRUE)
  df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)

  if(use.unadjusted.pvals) {
    ods <- which(df$lp<log(alpha))
  } else {
    ods <- which(df$lpa<log(alpha))
  }
  if(verbose) {
    print(paste0(length(ods), ' overdispersed genes ... ' ))
  }

  df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
  df$gsf[!is.finite(df$gsf)] <- 0;

  if(plot) {
    if(do.par) {
      par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
    }
    smoothScatter(df$m,df$v,main='',xlab='log10[ magnitude ]',ylab='log10[ variance ]')
    grid <- seq(min(df$m[vi]),max(df$m[vi]),length.out=1000)
    lines(grid,predict(m,newdata=data.frame(m=grid)),col="blue")
    if(length(ods)>0) {
      points(df$m[ods],df$v[ods],pch='.',col=2,cex=1)
    }
    smoothScatter(df$m[vi],df$qv[vi],xlab='log10[ magnitude ]',ylab='',main='adjusted')
    abline(h=1,lty=2,col=8)
    if(is.finite(max.adjusted.variance)) { abline(h=max.adjusted.variance,lty=2,col=1) }
    points(df$m[ods],df$qv[ods],col=2,pch='.')
  }

  ## variance normalize
  norm.mat <- counts*df$gsf
  if(!details) {
    return(rownames(mat)[ods])
  } else {
    ## return normalization factor
    return(list(mat=norm.mat, ods=rownames(mat)[ods], df=df))
  }
}
