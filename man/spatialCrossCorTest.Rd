% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlationStatistics.R
\name{spatialCrossCorTest}
\alias{spatialCrossCorTest}
\title{Tests for significance of spatial cross correlation for two features using random label null model}
\usage{
spatialCrossCorTest(x, y, w, n = 1000, ncores = 1, plot = FALSE, ...)
}
\arguments{
\item{x}{Feature 1 value}

\item{y}{Feature 2 value}

\item{w}{Binary weight matrix}

\item{n}{Permutation iterations}

\item{ncores}{Number of cores for parallel processing}

\item{plot}{Plot permutated distribution}

\item{...}{Additional parameters to pass to histogram plotting}
}
\value{
Two-sided test p-value
}
\description{
Tests for significance of spatial cross correlation for two features using random label null model
}
\examples{
\dontrun{
data(mOB)
pos <- mOB$pos
w <- getSpatialNeighbors(pos)
gexp <- normalizeCounts(mOB$counts, log=FALSE, verbose=FALSE)
pv1 <- spatialCrossCorTest(gexp['Gpsm1',], gexp['Nrgn',], w)
pv2 <- spatialCrossCorTest(gexp['Gpsm1',], gexp['Glul',], w)
pv3 <- spatialCrossCorTest(gexp[1,], gexp[2,], w)
}
}
