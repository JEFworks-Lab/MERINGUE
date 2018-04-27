# Simulation: https://gis.stackexchange.com/questions/161887/significance-test-for-morans-i-using-monte-carlo-simulation
#
# Moran's I
# https://en.wikipedia.org/wiki/Moran's_I
#
# Data locations.
#
x <- 1:8 # x-coordinates of points
y <- 1:5 # y-coordinates of points
#
# Weights matrix.
#
ind <- expand.grid(i=1:length(x), j=1:length(y))
f <- function(i, j) {
  u <- min(3, sum(abs(ind[i, ] - ind[j, ])))
  c(0, 1, sqrt(1/2), 0)[u+1]
}
w <- matrix(0.0, nrow(ind), nrow(ind))
for (i in 1:nrow(ind)) for (j in 1:nrow(ind)) w[i,j] <- f(i,j)
w <- w / sum(w)

#
# Simulated observations.
#
set.seed(17)
par(mfrow=c(2,3))

# Induce autocorrelation via an underlying trend.
z <- matrix(rexp(length(x)*length(y), outer(x,y^2)), length(x))
image(log(z), main="Autocorrelated Data")
hist(z)
library(reshape2)
zm <- melt(z)
pos <- zm[,1:2]
zv <- zm[,3]
plot(pos, col=map2col(zv), pch=16)

# double check that matrix and vector representation give same results
moranPermutationTest(zv, w, main="Null Distribution of Moran's I", xlab="I")
moranPermutationTest(z, w, main="Null Distribution of Moran's I", xlab="I")
# compare with empirical
moranTest(z, w)
moranTest(zv, w)

# Generate data independently of location.
z <- matrix(rnorm(length(x)*length(y), 0, 1/2), length(x))

