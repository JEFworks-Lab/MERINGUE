% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotNetwork}
\alias{plotNetwork}
\title{Plot an adjacency weight matrix as a network
Adapted from https://stackoverflow.com/questions/43879347/plotting-a-adjacency-matrix-using-pure-r}
\usage{
plotNetwork(pos, adj, col = "black", line.col = "grey", line.power = 1, ...)
}
\arguments{
\item{pos}{Position matrix}

\item{adj}{Adjacency weight matrix}

\item{col}{Color of points}

\item{line.col}{Color of line}

\item{line.power}{Thickness of lines}

\item{...}{Additional plotting parameters}
}
\value{
None
}
\description{
Plot an adjacency weight matrix as a network
Adapted from https://stackoverflow.com/questions/43879347/plotting-a-adjacency-matrix-using-pure-r
}
\examples{
data(mOB)
pos <- mOB$pos
w <- getSpatialNeighbors(pos)
plotNetwork(pos, w)

}
