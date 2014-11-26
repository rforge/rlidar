#' LiDAR - getVolume 3D
#'
#'@description This function calculates the volume of the ??-shape of the LiDAR point cloud in the three-dimensional space using alphashape3d package.
#'
#'@usage getVolume(xyz,id,alpha)
#'
#'@param xyz A 3-column matrix with the x,y and z coordinates of the LiDAR input points
#'@param alpha A value of the LiDAR input point id
#'@param alpha A single value or vector of values for ??.
#'@param cas Logical, if TRUE (default) plot the ??-shape.
#'@return Return dataframe of the LAS data set
#'@author Carlos Alberto Silva. Uses code by Beatriz Pateiro-Lopez (alphashape3d package)
#'@references \link{http://cran.r-project.org/web/packages/alphashape3d/index.html}
#'@examples
#'
#'
#'\dontrun{
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
#'#' set the cas parameter
#'cas=TRUE
#'
#'#' get the volume of the ??-shape
#'volume<-getVolume(xyz=xyz,id=id,cas=cas)
#'
#'}
#' @export
getVolume<-function(xyz,id,alpha,cas=TRUE) {

if (nrow(xyz)!=length(id)) {stop("The xyz and id do not have the same length")}

  xrange1<-range(xyz[,1])
  yrange2<-range(xyz[,2])
  ntree<-as.factor(id)
  N<-nlevels(ntree)
  
  xyzid<-cbind(xyz,id)
  
  ########FUNCTION TO RESCALED THE DATA############
  ## matlab equivalent of repmat in R
  repmat = function(X,m,n){
    ## tiles matrix X, m by n times
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
  
  # Scale to [0,1]x[0,1]x[0,1]
  # ==============================================
  real_range=diff(apply(allc,2,range))
  all_scaled1=allc-repmat(t(as.matrix(apply(allc,2,min))),nrow(allc),1)
  all_scaled=all_scaled1/repmat(t(as.matrix(apply(all_scaled1,2,max))),nrow(all_scaled1),1)
  
  # Compile this chunk for the scaled to [0,1]x[0,1]x[0,1] version
  # ================================================================
  # Original points (complete) scaled
  alpha=0.25
  alpha <- rep(alpha,N) # Alpha (the same for each clip)
  cuts<-c(0,cumsum(len)) # cuts+1 indicates the rows where each clip starts
  v<-numeric()
  volume_final<-numeric()

  for (i in 1:N){
    obs<-all_scaled[(cuts[i]+1):cuts[i+1],]  # Clip_i (from the complete scaled dataset!!!)
    b=ashape3d(as.matrix(obs),alpha=alpha[i],pert = TRUE, eps = 1e-09) # Ashape of clip_i (scaled)
    borig<-b # Ashape of clip_i (original scale)
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
    volume3d.trees[i]<-(volume_final[i]) # Total volume in the priginal scale
  }
  volume3d.trees<-cbind(id=as.numeric(levels(factor(id))),V=round(volume3d.trees, digits=2))
  colnames(volume3d.trees)<-c("id","V")
  
 return(volume3d.trees)
  
}

