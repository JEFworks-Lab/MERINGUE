#' Spatial transcriptomics transcriptomics of the mouse olfactory bulb
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"mOB"

#' MERFISH data of the mouse pre-optic region for a female naive animal (FN7)
#'
#' @format List where 'mat' is a sparse matrix with columns as cells and rows as genes
#'                          where expression values have already been normalized by volume
#'                    'pos' is a data frame of x, y, z position values per cell
#'                          and brain position as 6 slice indices from anterior to posterior
#'
#' @source \url{https://science.sciencemag.org/content/362/6416/eaau5324/}
"mPOA"

#' Spatial transcriptomics transcriptomics of 4 breast cancer biopsy sections
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'                          and slice index for 4 consecutive slices
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"BCL"

#' Drosophila embryo aligned ISH from the in situ database (BDTNP).
#'
#' @format List where 'mat' is a matrix with columns as cells and rows as genes and
#'                    'pos' is a data frame of x, y, and z position values per cells
#'
#' @source \url{https://shiny.mdc-berlin.de/DVEX/}
"drosophila"

#' SlideSeq data of the Purkinje layer of the mouse cerebellum for one puck (Puck_180819_12)
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'                          and slice index for 4 consecutive slices
#'
#' @source \url{https://science.sciencemag.org/content/363/6434/1463}
"purkinje"

#' Receptor ligand list
#'
#' @format Data frame corresponding to ncomms8866-s3.xlsx in Ramilowski et al (Nature Communications 2015)
#'
#' @source \url{https://www.nature.com/articles/ncomms8866}

#'#' Visium 10X Spatial Transcriptomics data of an adult mouse brain coronal section (P56)
#'
#' @format List where 'filteredGenes' is a sparse matrix with columns as genes and rows as spot IDs
#'                          where expression values have already been normalized and genes filtered
#'                          to 1263 genes whose expression variance across the 2702 spots is higher
#'                          than transcriptome-wide expectations. mt genes also removed.
#'                    'tissueSpotRotation' is a data frame of x, y position values per spot
#'                          where the spot barcodes are the 2702 that overlap the tissue.
#'                          Coordinated have been adjusted to match the high resolution PNG image
#'                          and rotated such that they visually match the tissue image when plotted.
#'                    'cluster6CorrMtx' is the cross correlation matrix of the significantly spatially
#'                          variable expressed genes in cluster 6.
#'
#' @source \url{https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/}
"mouseCoronal"

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
