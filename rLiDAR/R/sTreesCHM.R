#'Individual tree detection on the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@description Getting the individual tree location and hieght from the the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@usage sTreesCHM(chm, fws,htd)
#'
#'@param chm A LiDAR-derived Canopy Height Model (CHM)
#'@param fws Single dimension of fixed square window size, e.g. 3, 5, 7 and so on. Default is 5. 
#'@param htd Detection individual tree above specified heightbreak, e.g. 1.37, 2.0, 3.5 m and so on. Default is 1.37 m.
#'@return returns Individual tree dection list (x,y,and height)
#'@author Carlos Alberto Silva
#'@examples
#'
#'
#'\dontrun{
#'
#'#' Importing the LiDAR-derived CHM file
#'data(chm) # or set a CHM. e.g. chm<-readGDAL("CHM_stand.asc") 
#'
#'#' Set the fws:
#'fws<-5 # dimention 3x3
#'
#'#' Set the specified heightbreak
#'htd<-8.0
#'
#'#' Getting the individual tree detection list
#'treeList<-sTreesCHM(chm, fws,htd)
#'summary(treeList)
#'
#'#' Plotting the individual tree location
#'plot(chm) # plotting CHM
#'plot(SpatialPoints(treeList[,1:2]), add=T, col="red") # plotthing tree location
#'}
#' 
#' @export
sTreesCHM<-function(chm, fws=5,htd=1.37) {
  
  w<-matrix(c(rep(1,fws*fws)),nrow=fws,ncol=fws)
  
  if (class(chm)[1]!='RasterLayer') {chm<-raster(chm)}
  
  chm[chm < htd]<-NA
  
  f <- function(chm) max(chm) # , na.rm=TRUE
  rlocalmax <- focal(chm, fun=f, w=w, pad=TRUE, padValue=NA)
  
  setNull<- chm==rlocalmax
  XYmax <- SpatialPoints(xyFromCell(setNull, Which(setNull==1, cells=TRUE)))
  
  projection(XYmax)<-projection(chm)
  htExtract<-over(XYmax,as(chm, "SpatialGridDataFrame"))
  treeList<-cbind(XYmax,htExtract)
  
  colnames(treeList)<-c("x","y","height")
    
  return(treeList)
}
