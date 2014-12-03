#'@export
singleTreeCHM<-function(chm, fws=5,htd=1.37) {
  
  w<-matrix(c(rep(1,fws*fws)),nrow=fws,ncol=fws)
  
  if (class(chm)[1]!='RasterLayer') {chm<-raster(chm)}
  
  chm[chm < htd]<-NA
  
  f <- function(chm) max(chm)
  rlocalmax <- focal(chm, fun=f, w=w, pad=TRUE, padValue=NA)
  
  setNull<- chm==rlocalmax
  XYmax <- SpatialPoints(xyFromCell(setNull, Which(setNull==1, cells=TRUE)))
  
  projection(XYmax)<-projection(chm)
  htExtract<-over(XYmax,as(chm, "SpatialGridDataFrame"))
  treeList<-cbind(XYmax,htExtract)
  
  colnames(treeList)<-c("x","y","height")
    
  return(treeList)
}
