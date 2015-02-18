#'LiDAR - getVolume 3D
#'
#'@description This function calculates the volume of the 3D \eqn{\alpha}-shape of the LiDAR point cloud.
#'
#'@usage getVolume(xyz, id, alpha, plotCAS)
#'
#'@param xyz A 3-column matrix or dataframe with the x, y and z coordinates of the 3D LiDAR point cloud.
#'@param id A vector id for the x, y, z 3-column matrix. 
#'@param alpha A single value or vector of values for \eqn{\alpha}. Its range from 0 to 1.
#'@param plotCAS Logical, if TRUE (default) plot the \eqn{\alpha}-shape 3D.
#'@return Return volume of the 3D \eqn{\alpha}-shape of the LiDAR point cloud in units of cubic meters.
#'@author Carlos Alberto Silva. Uses code by Beatriz Pateiro-Lopez (\emph{alphashape3d} package,see \code{\link[alphashape3d]{volume_ashape3d}})
#'@references Lafarge, T.; Pateiro-Lopez, B.; Possolo, A. and Dunkers, J.P. (2014). R Implementation of a Polyhedral Approximation to a 3D Set of Points Using the alpha-Shape Journal of Statistical Software Vol. 56(4), pp. 1-18.
#'@examples
#'
#'\dontrun{
#'
#'# Importing LAS file:
#'LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")
#'
#'# Reading LAS file
#'LAS<-readLAS(LASfile,short=TRUE)
#'
#'# Setring the xyz coordinates and subsetting the data
#'xyz<-subset(LAS[,1:3],LAS[,3] >= 1.37)
#'
#'# Finding clusters
#'clLAS<-kmeans(xyz, 32)
#'
#'# set the id vector
#'id<-as.factor(clLAS$cluster)
#'
#'# set the alpha
#'alpha<-0.25
#'
#'# set the plotCAS parameter
#'plotCAS=TRUE
#'
#'# get the volume 
#'volume<-getVolume(xyz=xyz, id=id, alpha=alpha, plotCAS=plotCAS)
#'
#'# Add other plot parameters
#'aspect3d(1,1,0.5)
#'axes3d(c("x-", "x-", "y-", "z"), col="gray") # axes
#'title3d(xlab = "X Coord", ylab = " Y Coord", zlab = "Height", col="red") # title
#'
#'}
#'#@importFrom alphashape3d ashape3d volume_ashape3d
#'@importFrom rgl bg3d plot3d 
#'@export
getVolume<-function(xyz,id,alpha,plotCAS) {

  if (ncol(xyz)!=3) {stop("The input xyz must to be a 3-column matrix or dataframe with the x,y and z coordinates of the 3D LiDAR point cloud")}
  if (nrow(xyz)!=length(id)) {stop("The xyz and id do not have the same length")}
  if ((alpha>=0 & alpha<=1) !="TRUE") {stop("The alpha parameter is invalid. Please, use a value from 0 to 1")}
  if (class(plotCAS)!="logical") {stop("The plotCAS parameter is invalid. It must to be a TRUE or FALSE logical statement")}
  

  xrange1<-range(xyz[,1])
  yrange2<-range(xyz[,2])
  ntree<-as.factor(id)
  N<-nlevels(as.factor(ntree))
  
  xyzid<-cbind(xyz,id)
  
  repmat = function(X,m,n){
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  
  allc<-numeric(3)
  len=numeric()
  suppressWarnings(
  for (i in 1:N) {
    len[i]=dim(subset(xyzid,id==i))
    allc<-rbind(allc,subset(xyzid, id==i))
    
  })
  allc<-allc[-1,1:3]
  allc<-na.omit(rbind(allc,allc,allc))
  

  real_range=diff(apply(allc,2,range))
  all_scaled1=allc-repmat(t(as.matrix(apply(allc,2,min))),nrow(allc),1)
  all_scaled=all_scaled1/repmat(t(as.matrix(apply(all_scaled1,2,max))),nrow(all_scaled1),1)

  alpha <- rep(alpha,N) 
  cuts<-c(0,cumsum(len)) 
  v<-numeric()
  volume_final<-numeric()

  for (i in 1:N){
    obs<-all_scaled[(cuts[i]+1):cuts[i+1],] 
    b=ashape3d(as.matrix(obs),alpha=alpha[i],pert = TRUE, eps = 1e-09) 
    borig<-b 
    borig$x<-borig$x*repmat(t(as.matrix(apply(all_scaled1,2,max))),len[i],1)
    borig$x<-borig$x+repmat(t(as.matrix(apply(allc,2,min))),len[i],1)
    if (plotCAS==TRUE) {
      bg3d("white")
      plot(borig,indexAlpha = "all",clear=FALSE,col=rep(i,N),transparency = 1)}
    v[i]=volume_ashape3d(b)
    volume_final[i]=prod(v[i],real_range)
  }
  volume3d.trees<-NULL
  
  for( i in 1:N) {
    volume3d.trees[i]<-(volume_final[i]) 
  }
  volume3d.trees<-cbind(id=as.numeric(levels(factor(id))),V=round(volume3d.trees, digits=2))
  colnames(volume3d.trees)<-c("id","Volume")
 
  print(volume3d.trees)
  
  return(volume3d.trees)
 
 }

