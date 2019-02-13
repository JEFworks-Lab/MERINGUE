##### Tricking trendsceek with nonhomogenous poisson cell density
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library('devtools')
#devtools::install_github('edsgard/trendsceek')
library(trendsceek)
library(MERingue)

# simulate non-homogenous density with patterned expression at different length scales
set.seed(10)
pos <- cbind(
  c(rep(1,20), rep(2,20), rep(4, 20), rep(8, 20), rep(16, 20), rep(32, 20), rep(64, 20), rep(128, 20), rep(256, 20), rep(512, 20)),
  c(rnorm(20), rnorm(20), rnorm(20), rnorm(20), rnorm(20), rnorm(20), rnorm(20), rnorm(20), rnorm(20), rnorm(20))
)
colnames(pos) <- c('x', 'y')
pos[,1]=pos[,1]/max(pos[,1])
pos[,2]=pos[,2]/max(pos[,2])
gexp <- c(rep(0, 20), rep(1,20), rep(0,20), rep(1, 20), rep(0, 20), rep(1, 20), rep(0, 20), rep(1, 20), rep(0, 20), rep(1, 20))
#gexp <- jitter(scale(gexp)[,1])
rownames(pos) <- names(gexp) <- paste0('cell', 1:length(gexp))
plot(pos)
par(mfrow=c(2,1))
hist(pos[,1], breaks=1000)
plotEmbedding(pos, col=gexp, cex=2)

# create pp object for trendsceek
pp = sim_pois(nrow(pos))
pp$n <- length(gexp)
pp$x=pos[,1]
pp$y=pos[,2]
pp$marks=data.frame(gexp)
pp$marks
trendstat_list = trendsceek_test(pp, nrand = 100, ncores = 1)
do.call(cbind, lapply(trendstat_list[["supstats"]], function(x) x$p.bh))

## MERingue
w <- voronoiAdjacency(pos, plot=TRUE)
plotNetwork(pos, w)
#w <- getSpatialWeights(pos, 3)
#plotNetwork(pos, w)
do.call(cbind, MERingue:::moranTest(gexp, w, alternative='two.sided'))

# change back to more uniform density
pos[,1] <- log(pos[,1])
pos[,1]=pos[,1] - min(pos[,1])
pos[,1]=pos[,1]/max(pos[,1])
plot(pos)
par(mfrow=c(2,1))
hist(pos[,1], breaks=100)
plotEmbedding(pos, col=gexp, cex=2)

# repeat to get significant results
pp = sim_pois(nrow(pos))
pp$n <- length(gexp)
pp$x=pos[,1]
pp$y=pos[,2]
pp$marks=data.frame(gexp)
pp$marks
trendstat_list = trendsceek_test(pp, nrand = 100, ncores = 1)
trendstat_list[["supstats"]]
do.call(cbind, lapply(trendstat_list[["supstats"]], function(x) x$p.bh))

par(mfrow=c(1,1))
w <- voronoiAdjacency(pos, plot = TRUE)
do.call(cbind, MERingue:::moranTest(gexp, w, alternative='two.sided'))






############# Supp Fig 1.
## Simulating failure examples

## A. Random but misinterpreted as clustered due to cell density
expPos <- function(y) {
  x <- y
  a <- x[x<0]
  b <- x[x>=0]
  a <- -(log10(abs(a))+1)
  b <- log10(b+1)
  y[x<0] <- a
  y[x>=0] <- b
  return(y)
}
set.seed(0)
k <- 100
gexp <- rnorm(k)
pos <- cbind(expPos(rnorm(k, 100, 100)), expPos(rnorm(k, 100, 100)))
#pos <- cbind(rnorm(k, 100, 100), rnorm(k, 100, 100))

#pos[,1]=pos[,1] - min(pos[,1])
#pos[,1]=pos[,1]/max(pos[,1])

#pos[,2]=pos[,2] - min(pos[,2])
#pos[,2]=pos[,2]/max(pos[,2])

range(pos)

rownames(pos) <- names(gexp) <- 1:k
MERingue::plotEmbedding(pos, colors=gexp, cex=1)


library(trendsceek)
pp = pos2pp(pos)
log.fcn = log10
pp = set_marks(pp, t(10^gexp), log.fcn = log.fcn)
nrand = 100
ncores = 1
trendstat_list = trendsceek_test(pp, nrand, ncores)
trendstat_list[["supstats"]]
do.call(cbind, lapply(trendstat_list[["supstats"]], function(x) x$p.bh))

par(mfrow=c(1,1))
w <- voronoiAdjacency(pos, plot = TRUE)
do.call(cbind, MERingue:::moranTest(gexp, w, alternative='two.sided'))














par(mfrow=c(1,1))

set.seed(0)
N <- 10^2
pos <- t(combn(c(1:sqrt(N), rev(1:sqrt(N))), 2))
pos <- unique(pos)
dim(pos)
rownames(pos) <- paste0('cell', 1:N)
colnames(pos) <- c('x', 'y')
dim(pos)
plotEmbedding(pos)

pos[,1] <- expPos(pos[,1])
pos[,2] <- expPos(pos[,2])
plot(pos)

gexp <- rep(0, N)
#vi <- pos[,1] < 0.5 | pos[,2] < 0.5
a <- unique(c(pos[,1], pos[,2]))
foobar <- a[seq(1, length(a), 2)]
foobar
#vi <- !(pos[,1] %in% foobar | pos[,2] %in% foobar)
vi <- pos[,1] %in% foobar
table(vi)
gexp[vi] <- 1

rownames(pos) <- names(gexp) <- paste0('cell', 1:length(gexp))
colnames(pos) <- c('x','y')

plotEmbedding(pos, colors=jitter(scale(gexp)[,1]), cex=1)


library(trendsceek)
pp = pos2pp(pos)
log.fcn = log10
pp = set_marks(pp, t(10^gexp), log.fcn = log.fcn)
trendstat_list = trendsceek_test(pp, 100, 1)
trendstat_list[["supstats"]]
do.call(cbind, lapply(trendstat_list[["supstats"]], function(x) x$p.bh))

par(mfrow=c(1,1))
w <- voronoiAdjacency(pos, plot = TRUE)
do.call(cbind, MERingue:::moranTest(gexp, w, alternative='two.sided'))
