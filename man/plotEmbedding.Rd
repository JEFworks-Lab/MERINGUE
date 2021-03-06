% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotEmbedding}
\alias{plotEmbedding}
\title{Plot 2D embedding}
\usage{
plotEmbedding(
  emb,
  groups = NULL,
  colors = NULL,
  cex = 0.6,
  alpha = 0.4,
  gradientPalette = NULL,
  zlim = NULL,
  s = 1,
  v = 0.8,
  min.group.size = 1,
  show.legend = FALSE,
  mark.clusters = FALSE,
  mark.cluster.cex = 2,
  shuffle.colors = FALSE,
  legend.x = "topright",
  gradient.range.quantile = 0.95,
  verbose = FALSE,
  unclassified.cell.color = "lightgrey",
  group.level.colors = NULL,
  xlab = NA,
  ylab = NA,
  ...
)
}
\arguments{
\item{emb}{dataframe with x and y coordinates}

\item{groups}{factor annotations for rows on emb for visualizing cluster annotations}

\item{colors}{color or numeric values for rows on emb for visualizing gene expression}

\item{cex}{point size}

\item{alpha}{point opacity}

\item{gradientPalette}{palette for colors if numeric values provided}

\item{zlim}{range for colors}

\item{s}{saturation of rainbow for group colors}

\item{v}{value of rainbow for group colors}

\item{min.group.size}{minimum size of group in order for group to be colored}

\item{show.legend}{whether to show legend}

\item{mark.clusters}{whether to mark clusters with name of cluster}

\item{mark.cluster.cex}{cluster marker point size}

\item{shuffle.colors}{whether to shuffle group colors}

\item{legend.x}{legend position ie. 'topright', 'topleft', 'bottomleft', 'bottomright'}

\item{gradient.range.quantile}{quantile for mapping colors to gradient palette}

\item{verbose}{verbosity}

\item{unclassified.cell.color}{cells not included in groups will be labeled in this color}

\item{group.level.colors}{set group level colors. Default uses rainbow.}

\item{xlab}{x-axis label.}

\item{ylab}{y-axis label.}

\item{...}{Additional parameters to pass to BASE::plot}
}
\value{
None
}
\description{
Plot 2D embedding
}
\examples{
data(mOB)
pos <- mOB$pos
gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)['Camk4',]
plotEmbedding(pos, colors=scale(gexp)[,1], zlim=c(-2,2), cex=3)

}
