#' LiDAR - getVolume 3D
#'
#'@description This function calculates the volume of the alpha-shape of the LiDAR point cloud in the three-dimensional space using alphashape3d package.
#'
#'@usage getVolume(xyz,id,alpha,cas)
#'
#'@param xyz A 3-column matrix with the x,y and z coordinates of the LiDAR input points
#'@param id A vector id for the the items in xyz. 
#'@param alpha A single value or vector of values for alpha. Its range from 0 to 1.
#'@param cas Logical, if TRUE (default) plot the alpha-shape.
#'@return Return dataframe of the LAS data set
#'@author Carlos Alberto Silva. Uses code by Beatriz Pateiro-Lopez (alphashape3d R-package)
#'@references \link{http://cran.r-project.org/web/packages/alphashape3d/index.html}
#'@examples
#'
#'\dontrun{
#'
#'#' Importing LAS file:
#'myLAS<-data(LASfile) # or set a LAS  file (myLAS<-"LASfile.las")
#'
#'#' Reading LAS file
#'#'LAS<-readLAS(myLAS,short=TRUE)
#'
#'#' Finding clusters
#'clLAS<-kmeans(LAS[,3], 32)
#'
#'# set the xyz coordenates
#'xyz<-LAS[,1:3]
#'
#'#' set the id vector
#'id<-as.factor(clLAS$cluster)
#'
#' set the alpha value
#'alpha<-0.25
#'
#'#' set the cas parameter
#'cas=TRUE
#'
#'#' get the volume of the alpha-shape
#'volume<-getVolume(xyz=xyz,id=id,alpha=alpha,cas=cas)
#'head(volume)
#'}
#'
#' @export
getVolume<-function(xyz,id,alpha,cas) {

if (nrow(xyz)!=length(id)) {stop("The xyz and id do not have the same length")}
if (cas!=TRUE & cas!=FALSE) {stop("The cas parameter is invalid. Please, use TRUE or FALSE")}

  xrange1<-range(xyz[,1])
  yrange2<-range(xyz[,2])
  ntree<-as.factor(id)
  N<-nlevels(ntree)
  
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

  alpha=0.25
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
    if (cas==TRUE) {
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
  colnames(volume3d.trees)<-c("id","V")
  
 return(volume3d.trees)
 }

