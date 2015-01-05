#'Individual tree detection on the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@description Detect and compute individual individual tree location and height on the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@usage singleTreeCHM(chm,fws,minht)
#'
#'@param chm A LiDAR-derived Canopy Height Model (CHM) raster  file
#'@param fws A single dimension of fixed square window size, e.g. 3, 5, 7 and so on. Default is 5 
#'@param minht Detect individual tree above specified heightbreak, e.g. 1.37, 2.0, 3.5 m and so on. Default is 1.37 m
#'@return returns A 3-column matrix of the location and height for the individual trees dectected
#'@author Carlos Alberto Silva
#'@examples
#'\dontrun{
#'
#'# Importing the LiDAR-derived CHM raster file
#'
#'library(raster)
#'data(chm) # or set a CHM. e.g. chm<-raster("CHM_stand.asc") 
#'
#'# Smoothing CHM
#'schm<-CHMsmoothing(chm, "mean", 5)
#'
#'# Plotting smothed CHM
#'plot(schm) 
#'
#'# Set the fws:
#'fws<-5 # dimention 5x5
#'
#'# Set the specified heightbreak
#'minht<-8.0
#'
#'# Getting the individual tree detection list
#'treeList<-singleTreeCHM(schm, fws, minht)
#'summary(treeList)
#'
#'# Plotting the individual tree location on the CHM
#'library(sp)
#'plot(chm) # plotting CHM
#'plot(SpatialPoints(treeList[,1:2]), add=T, col="red") # plotthing tree location
#'}
#'
#'@importFrom raster raster focal xyFromCell Which projection
#'@importFrom sp SpatialPoints over SpatialGridDataFrame 
#'@export
singleTreeCHM<-function(chm, fws=5,minht=1.37) {
  
  if (class(chm)[1]!='RasterLayer') {chm<-raster(chm)}
  if (class(fws)!="numeric") {stop("The fws parameter is invalid. It is not a numeric input")}
  if (class(minht)!="numeric") {stop("The minht parameter is invalid. It is not a numeric input")}
  
  w<-matrix(c(rep(1,fws*fws)),nrow=fws,ncol=fws)
    
  chm[chm < minht]<-NA
  
  f <- function(chm) max(chm)
  
  rlocalmax <- focal(chm, fun=f, w=w, pad=TRUE, padValue=NA)
  
  setNull<- chm==rlocalmax
  XYmax <- SpatialPoints(xyFromCell(setNull, Which(setNull==1, cells=TRUE)))
  
  projection(XYmax)<-projection(chm)
  htExtract<-over(XYmax,as(chm, "SpatialGridDataFrame"))
  treeList<-cbind(XYmax,htExtract)
  
  colnames(treeList)<-c("x","y","height")
    
  return(treeList)
  print(treeList)
}
