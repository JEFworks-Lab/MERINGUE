# Code to support Toroidal shift model (TSM) as null model

# Code mofied from http://www.is.titech.ac.jp/~mase/mase/splancs/

# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#

bbox <- function(pts)
{

  xr <- range(pts[,1],na.rm=T)
  yr <- range(pts[,2],na.rm=T)

  cbind(c(xr[1],xr[2],xr[2],xr[1]),c(yr[1],yr[1],yr[2],yr[2]))
}

shift <- function(pts,xsh=0.0,ysh=0.0)
{
  pts[,1] <- pts[,1]+xsh
  pts[,2] <- pts[,2]+ysh
  pts
}

torShift <- function(pts,xsh,ysh,rect=NULL)
{
  if(is.null(rect)) { rect <- bbox(pts) }

  xoff <- min(rect[,1])
  yoff <- min(rect[,2])

  xsc <- (max(rect[,1])-xoff)
  ysc <- (max(rect[,2])-yoff)

  pts[,1] <- pts[,1]-xoff
  pts[,2] <- pts[,2]-yoff
  pts <- shift(pts,xsh,ysh)
  pts[,1] <- (pts[,1] %% xsc )+xoff
  pts[,2] <- (pts[,2] %% ysc )+yoff
  pts
}


rtorShift <- function(pts, rect=NULL, k=4, seed=0)
{
  set.seed(seed)
  if(is.null(rect)) { rect <- bbox(pts)/k }
  xsc <- max(rect[,1])-min(rect[,1])
  ysc <- max(rect[,2])-min(rect[,2])
  xsh <- runif(1)*xsc
  ysh <- runif(1)*ysc
  torShift(pts,xsh,ysh,rect)
}
