#'LiDAR individual tree canopy area 
#'
#'@description Compute and export individual tree canopy area detected on the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@usage ForestCAS(chm,loc,maxcrown,exclusion)
#'
#'@param chm A raster LiDAR-derived Canopy Height Model (CHM)
#'@param loc A 3-column matrix with the x,y coordinates and heights of the individual tree
#'@param maxcrown A single value of the maximum individual tree crown radius expected. Default 10.0 m
#'@param exclusion A single value with the percent of pixel exclusion. E.g 0.5 exclusion of the pixels that has values below of the 50 of max tree height. Default is 0.3 
#'@return returns A list that contains the individual tree canopy boundary polygon and the dataframe of the canopy area  
#'@author Carlos Alberto Silva
#'@examples
#'\dontrun{
#'# Importing the LiDAR-derived CHM file
#'data(chm) # or set a CHM. e.g. chm<-raster("CHM_stand.asc") 
#'
#'# Set the loc parameter
#'sCHM<-CHMsmoothing(chm, filter="mean", ws=5, sigma=NULL) # smoothing CHM
#'loc<-singleTreesCHM(sCHM, fws=5,htd=8) # or import a tree list
#'
#'# Set the maxcrown parameter
#'maxcrown=10.0 
#'
#'# Set the exclusion parameter
#'exclusion=0.3 # 30
#'
#'# Compute individual tree detection canopy area
#'canopy<-ForestCAS(chm,loc,maxcrown,exclusion)
#'
#'# Getting the individual tree detection canopy area boundary
#'boundaryTrees<-canopy[[1]]
#
#'# Plotting the individual tree canopy boundary over the CHM
#'plot(chm) # plotting CHM
#'plot(boundaryTrees, add=T, border='red', bg='transparent') # adding tree canopy boundary
#'
#'# Getting the individual tree detection canopy area list
#'canopyList<-canopy[[2]]
#'summary(canopyList)
#'plot(SpatialPoints(canopyList[,1:2]),col="black", add=T, pch="*") # adding tree location to the plot
#'} 
#'@export
ForestCAS<-function(chm,loc,maxcrown,exclusion) {

  chm<-as(chm, "SpatialGridDataFrame")
  Hthreshold<-min(loc[,3])*exclusion
  polys<-list() 
  width<-numeric()  
  
  for(i in 1:nrow(loc)) { 
    width[i] =maxcrown
    discbuff<-disc(radius=width[i], centre=c(loc$x[i], loc$y[i])) 
    discpoly<-Polygon(rbind(cbind(discbuff$bdry[[1]]$x, 
                                  y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1], 
                                                             y=discbuff$bdry[[1]]$y[1]))) 
    polys<-c(polys, discpoly) 
  } 
  
  spolys<-list() 
  for(i in 1:length(polys)) { 
    spolybuff<-Polygons(list(polys[[i]]), ID=row.names(loc)[i]) 
    spolys<-c(spolys, spolybuff) 
  } 
  polybuffs<-SpatialPolygons(spolys) 

  chmdf<-as.data.frame(chm,xy=TRUE)
  Points.Ply<-over(SpatialPoints(chmdf[,2:3]),polybuffs) # overlay
  Points.PlyD<-cbind(chmdf,Points.Ply) # cbind
  Points.PlyD<-na.omit(Points.PlyD) # omite NAs
 vor =  deldir(loc[,1], loc[,2], z=loc[,3],suppressMsge=T)
  tile = tile.list(vor)
  polys = vector(mode='list', length=length(tile))
  for (i in seq(along=polys)) {
    pcrds = cbind(tile[[i]]$x, tile[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  
  SP = SpatialPolygons(polys)
  veronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=loc[,1], 
                                                         y=loc[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                          function(x) slot(x, 'ID'))))
  chmV<-over(SpatialPoints(Points.PlyD[,2:3]), SP)
  RpD<-cbind(Points.PlyD,chmV)
  RpD.filter<-subset(RpD[,1:5],RpD[,1]>=Hthreshold)
  RpD.filter<-cbind(RpD.filter[,1:3],RpD.filter[,5])
  colnames(RpD.filter)<-c("z","x","y","g")
  h.mH<-ddply(RpD.filter,.(g), function (RpD.filter)
    subset(RpD.filter,RpD.filter[,1]>= max(RpD.filter[,1])*exclusion))
  
  DF2raster<-function(h.mH, i){
    h.mHkl<-subset(h.mH,h.mH[,4]==levels(factor(h.mH[,4]))[i])
    
    if (nrow(h.mHkl==1)==TRUE) {
      h.mHkl<-rbind(h.mHkl,c(h.mHkl[,1],h.mHkl[,2]+0.005,h.mHkl[,3]+0.005,h.mHkl[,4]))}
    
    spP <- cbind(h.mHkl[,2:3],h.mHkl[,1],h.mHkl[,4])
    coordinates(spP)<- ~ x + y
    suppressWarnings(gridded(spP) <- TRUE)
    
    rasterDF <- raster(spP)
    hhg<- boundaries(rasterDF, type='outer') 
    p <- rasterToPolygons(hhg, dissolve=TRUE)
    sp.polys <- p[1,]
    return(sp.polys)}
 
  for ( j in 1:nlevels(factor(h.mH[,4]))){
    assign(paste0("SP.polys", j), DF2raster(h.mH,j))
    print(paste("computting canopy area: Tree",j))}
  
  polygons <- slot(SP.polys1, "polygons")
  
  for (i in 1:nlevels(factor(h.mH[,4]))) {
    data.loc <- get(paste0("SP.polys",i))
    polygons <- c(slot(data.loc, "polygons"),polygons)
  }
  
  for (i in 1:length(polygons)) {
    slot(polygons[[i]], "ID") <- paste(i)
  }
  
  spatialPolygons <- SpatialPolygons(polygons)
  spdf <- SpatialPolygonsDataFrame(spatialPolygons, 
                                   data.frame(Trees=1:length(polygons)))
  
  options(scipen=10)
  spdf<-spdf[spdf@data[-length(polygons),],]
  areaList<-sapply(slot(spdf, "polygons"), slot, "area")
  canopyTable<-cbind(loc,areaList)
  colnames(canopyTable)<-c("x","y","z","ca")
  result=list(spdf,canopyTable) 
  return(result)
  
}