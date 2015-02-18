DF2raster<-function(h.mH, i){
  
  h.mHkl<-subset(h.mH,h.mH[,4]==levels(factor(h.mH[,4]))[i])
  
  if (nrow(h.mHkl==1)==TRUE) {
    h.mHkl<-rbind(h.mHkl,c(h.mHkl[,1],h.mHkl[,2]+0.005,h.mHkl[,3]+0.005,h.mHkl[,4]))}
  
  spP <- cbind(h.mHkl[,2:3],h.mHkl[,1],h.mHkl[,4])
  colnames(spP)<-c("x","y","z","y")
  #coordinates(spP)<- c("x", "y")
  #suppressWarnings(gridded(spP) <- TRUE)
  m = suppressWarnings(SpatialPixelsDataFrame(points=spP[c("x", "y")], data = spP))
  
  rasterDF <- raster(m)
  hhg<- boundaries(rasterDF, type='outer') 
  p <- rasterToPolygons(hhg, dissolve=TRUE)
  sp.polys <- p[1,]
  return(sp.polys)
}