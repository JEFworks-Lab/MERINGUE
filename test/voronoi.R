# Reanalysis
cd <- read.csv('../../../SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('../../../SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]
pos <- pos[,1:2]

# Voronoi
#install.packages('deldir')
library(deldir)
vtess <- deldir(pos[,1], pos[,2])

par(mfrow=c(1,1))
plot(pos)
plot(vtess, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1)

vor_desc <- tile.list(vtess) # tile.list extracts the polygon data from the deldir computation
# ^^ gets us the points for the polygons but we still have to close them, hence the need for the rbind

library(sp)
vor_polygons <- lapply(1:(length(vor_desc)), function(i) {

    tmp <- cbind(vor_desc[[i]]$x, vor_desc[[i]]$y)
    tmp <- rbind(tmp, tmp[1,])

    Polygons(list(Polygon(tmp)), ID=i) # now we can make the Polygon(s)

})

# dummy data frame that makes IDs easier
xdf <- data.frame(id=sapply(slot(SpatialPolygons(vor_polygons), 'polygons'), slot, 'ID'))
rownames(xdf) <- xdf$id

p <- SpatialPolygonsDataFrame(SpatialPolygons(vor_polygons), data=xdf)

library(spdep)
xy <- pos
w <- poly2nb(p, row.names=p$Id)
plot(p, col='gray', border='blue', lwd=2)
plot(w, xy, col='red', lwd=2, add=TRUE)

ww <-  nb2listw(w, style='B')
source('../R/process.R')
mat <- normalizeCounts(t(as.matrix(cd)))
gexp <- mat[1,]
moran(gexp, ww, n=length(ww$neighbours), S0=Szero(ww))
