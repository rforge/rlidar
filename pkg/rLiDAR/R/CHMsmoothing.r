#'LiDAR-derived Canopy Height Model (CHM) smoothing
#'
#'@description LiDAR-derived Canopy Height Model (CHM) smooth using focal statistcs
#'
#'@usage CHMsmoothing(chm, filter, ws, sigma)
#'
#'@param chm A LiDAR-derived Canopy Height Model (CHM) raster file.
#'@param ws A single square matrix of weights dimension, e.g. 3,5, 7 and so on. Default is 5.
#'@param filter mean, median, maximum, minimum or gaussian. Default is "mean".
#'@param sigma Used only when filter parameter is equal to gaussian, e.g. 0.5, 1.0, 1.5 and so on. Default is 0.6. 
#'@return returns A CHM smoothed raster.
#'@author Carlos Alberto Silva. 
#'@seealso Focal in R raster package at \code{\url{http://cran.r-project.org/web/packages/raster/}}.
#'@examples
#'\dontrun{
#'
#'# Importing the LiDAR-derived CHM file
#'data(chm) # or set a CHM. e.g. chm<-raster("CHM_stand.asc") 
#'
#'#------------------------------------------------#
#'# Example 01: Smoothing CHM using gaussian filter
#'
#'# Set the ws:
#'ws<-3 # dimention 3x3
#'
#'# Set the filter type
#'filter<-"gaussian"
#'
#'# Set the sigma value
#'sigma<-0.6
#'
#'# Smoothing CHM
#'sCHM<-CHMsmoothing(chm, filter, ws, sigma)
#'
#'# Plotting CHM smoothed
#'plot(sCHM, main=paste(filter,"filter and windows size", paste0(ws,"x",ws)))
#'
#'#------------------------------------------------# 
#'# Example 02: Smoothing CHM using mean filter
#'
#'# Set the ws:
#'ws<-5 # dimention 5x5
#'
#'# Set the filter type
#'filter<-"mean"
#'
#'# Smoothing CHM
#'sCHM<-CHMsmoothing(chm, filter, ws, sigma=NULL)
#'
#'# Plotting CHM smoothed
#'plot(sCHM, main=paste(filter,"filter and window size", paste0(ws,"x",ws)))
#'}
#'
#'@importFrom raster raster focal
#'@export CHMsmoothing
CHMsmoothing<-function(chm, filter="mean", ws=5, sigma=0.6) {

  if (class(chm)[1]!='RasterLayer') {
      chmInput<-raster(chm)
      } else {chmInput<-chm
  }
  
  if (filter == "mean") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=mean)
  }
  if (filter == "median") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=median)
  }
  if (filter == "maximum") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=max)
  }
  if (filter == "minimum") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=min)
  }
  
  if (filter =="gaussian") {
    
    fgauss <- function(sigma, n=ws) {
      m <- matrix(ncol=n, nrow=n)
      col <- rep(1:n, n)
      row <- rep(1:n, each=n)
      x <- col - ceiling(n/2)
      y <- row - ceiling(n/2)
      m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
      m / sum(m)
    }
    gf=fgauss(sigma)
    chmR <- focal(chmInput, w=gf)
  }
  return(chmR)
}
