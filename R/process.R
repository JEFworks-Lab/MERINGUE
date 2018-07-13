#' Filter a counts matrix
#'
#' @description Filter a counts matrix based on gene (row) and cell (column)
#'      requirements.
#'
#' @param counts A sparse read count matrix. The rows correspond to genes,
#'      columns correspond to individual cells
#' @param min.lib.size Minimum number of genes detected in a cell. Cells with
#'      fewer genes will be removed (default: 1)
#' @param max.lib.size Maximum number of genes detected in a cell. Cells with
#'      more genes will be removed (default: Inf)
#' @param min.reads Minimum number of reads per gene. Genes with fewer reads
#'      will be removed (default: 1)
#' @param min.detected Minimum number of cells a gene must be seen in. Genes
#'      not seen in a sufficient number of cells will be removed (default: 1)
#' @param verbose Verbosity (default: TRUE)
#'
#' @return a filtered read count matrix
#'
#' @export
#'
#' @importFrom Matrix Matrix colSums rowSums
#'
cleanCounts <- function (counts, min.lib.size = 1, max.lib.size = Inf, min.reads = 1, min.detected = 1, verbose = TRUE) {
  if (!class(counts) %in% c("dgCMatrix", "dgTMatrix", "dgeMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (verbose) {
    message("Filtering matrix with ", ncol(counts), " cells and ",
            nrow(counts), " genes ...")
  }
  ix_col <- Matrix::colSums(counts)
  ix_col <- ix_col > min.lib.size & ix_col < max.lib.size
  counts <- counts[, ix_col]
  counts <- counts[Matrix::rowSums(counts) > min.reads, ]
  counts <- counts[Matrix::rowSums(counts > 0) > min.detected, ]
  if (verbose) {
    message("Resulting matrix has ", ncol(counts), " cells and ", nrow(counts), " genes")
  }
  return(counts)
}


#' Normalizes counts to CPM
#'
#' @description Normalizes raw counts to log10 counts per million with pseudocount
#'
#' @param counts Read count matrix. The rows correspond to genes, columns
#'      correspond to individual cells
#' @param normFactor Normalization factor such as cell size. If not provided
#'      column sum as proxy for library size will be used
#' @param depthScale Depth scaling. Using a million for CPM (default: 1e6)
#' @param pseudo Pseudocount for log transform (default: 1)
#' @param verbose Verbosity (default: TRUE)
#'
#' @return a normalized matrix
#'
#' @export
#'
#' @importFrom Matrix Matrix colSums t
#'
normalizeCounts <- function (counts, normFactor = NULL, depthScale = 1e+06, pseudo=1, verbose = TRUE) {
  if (!class(counts) %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (verbose) {
    message("Normalizing matrix with ", ncol(counts), " cells and ", nrow(counts), " genes.")
  }
  if(is.null(normFactor)) {
    if (verbose) {
      message('normFactor not provided. Normalizing by library size.')
    }
    normFactor <- Matrix::colSums(counts)
  }
  if (verbose) {
    message(paste0("Using depthScale ", depthScale, " and log transforming with pseudocount ", pseudo,"."))
  }
  counts <- Matrix::t(Matrix::t(counts)/normFactor)
  counts <- counts * depthScale
  counts <- log10(counts + pseudo)
  return(counts)
}
