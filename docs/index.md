![](tools/img/meringue_logo.svg)
# MERINGUE

`MERINGUE` enables spatial gene expression analysis in non-homogenous tissues. The overall approach is detailed in the following publication: **COMING SOON**

## Benefits and Capabilities

(1) Provide a statistical framework to identify and characterize significantly spatially variable genes

(2) Group significantly spatially variable genes into primary spatial gene expression patterns

(3) Test for putative cellular interactions between spatially co-localized cell-types

(4) Perform spatially-informed transcriptional clustering to identify spatially-distinct cell-types and cell-states

(5) Accomodates 2D, multi-section, and 3D spatial data 

(6) Is robut to variations in cellular densities, distortions, or warping common to tissues

(7) Highly scalable to enable analysis of 10,000s of genes and 1,000s of cells within seconds

## Installation

To install `MERINGUE`, we recommend using `devtools`:
```
require(devtools)
devtools::install_github('JEFworks/MERINGUE')
```
## Tutorials

[mOB Spatial Transcriptomics Analysis](mOB_analysis)

[Multi-section 3D Breast Cancer Spatial Transcriptomics Analysis](BCL_analysis)

[3D Drosophila Spatial Transcriptomics Analysis](drosophila_3D_analysis)

[Understanding MERINGUE's Spatial Cross-Correlation Statistic using Simulations](simulation)

[Spatially-informed transcriptional clustering with MERINGUE](spatiall_clustering)

## Contributing

We welcome any bug reports, enhancement requests, general questions, and other contributions. To submit a bug report or enhancement request, please use the `MERINGUE` GitHub issues tracker. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.
