#' R6 object for storing spatial expression data
#'
SpatialExpressionData <- R6::R6Class(

  "SpatialExpressionData",

  public = list(
    name = NULL,
    exprs = NULL,
    pos = NULL,
    annot = NULL,
    verbose = NULL,

    initialize = function(name = NA, exprs = NA, pos = NA, annot = NULL, verbose = TRUE, ...) {
        self$name <- name
        self$exprs <- exprs

        vi <- intersect(rownames(pos), colnames(exprs))
        if(length(vi) < ncol(exprs)) {
          warning(paste0(length(vi), ' expression datasets with positional information'))
        }
        self$pos <- pos[vi,]

        self$annot <- annot
        self$verbose <- verbose
    }
  )

)

