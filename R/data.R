#' Mouse olfactory bulb spatial transcriptomics transcriptomics
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"mOB"

#' Breast cancer biopsy layer 1 bulb spatial transcriptomics transcriptomics
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"BCL1"

#' Breast cancer biopsy layer 2 bulb spatial transcriptomics transcriptomics
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"BCL2"

#' Breast cancer biopsy layer 3 bulb spatial transcriptomics transcriptomics
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"BCL3"

#' Breast cancer biopsy layer 4 bulb spatial transcriptomics transcriptomics
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"BCL4"

#' Drosophila embryo aligned ISH
"drosophila"

#' SlideSeq data
"purkinje"

#' Receptor ligand list
#'
#' @format Data frame corresponding to ncomms8866-s3.xlsx in Ramilowski et al (Nature Communications 2015)
#'
#' @source \url{https://www.nature.com/articles/ncomms8866}
#'
#' @examples
#' data(receptorLigandInfo)
#' receptors <- unique(receptorLigandInfo$Receptor.ApprovedSymbol)
#' receptorLigandList <- lapply(receptors, function(rp) {
#'      vi <- receptorLigandInfo$Receptor.ApprovedSymbol == rp
#'      return(receptorLigandInfo[vi,]$Ligand.ApprovedSymbol)
#' })
#' names(receptorLigandList) <- receptors
#' ligands <- unique(receptorLigandInfo$Ligand.ApprovedSymbol)
#' ligandReceptorList <- lapply(ligands, function(rp) {
#'      vi <- receptorLigandInfo$Ligand.ApprovedSymbol == rp
#'      return(receptorLigandInfo[vi,]$Receptor.ApprovedSymbol)
#' })
#' names(ligandReceptorList) <- ligands
"receptorLigandInfo"
