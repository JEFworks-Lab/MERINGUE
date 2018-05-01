# Reanalysis
cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# Voronoi
install.packages('deldir')
library(deldir)
vtess <- deldir(pos[,1], pos[,2])

par(mfrow=c(1,1))
plot(pos)
plot(vtess, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1)

