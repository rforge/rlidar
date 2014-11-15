#'LiDAR 3D stand visualization of trees
#'
#'@description LiDAR 3D stand visualization of trees. This function was adapted from MAESPA package,see
#'\code{\link{MAESPA}}.  
#'
#'@usage LiDARtrees3D(crownshape = c("cone", "ellipsoid", "halfellipsoid",
#'                 "paraboloid", "cylinder"), CL = 4, CW = 8, HCB = 10, 
#'                   X = 0, Y = 0, dbh = 0.3, crowncolor = "forestgreen", 
#'                  stemcolor = "chocolate4", resolution="high", shape=TRUE) 
#'
#'@param crownshape Tree crown shape: "cone", "ellipsoid","halfellipsoid", "paraboloid" or "cylinder". Defaul is "halfellipsoid"
#'@param CL Tree crown length
#'@param CW Tree crown diameter
#'@param HCB Tree trunk height
#'@param X X Tree location
#'@param Y Y Tree location
#'@param dbh Tree stem diameter
#'@param crowncolor Tree crown color
#'@param stemcolor Tree stem color
#'@param resolution Tree crown resolution: "low", "median" and "high"
#'@param shape TRUE return a interpolate tree crown shape  
#'
#'@return returns Single 3D tree 
#'@author Carlos Alberto Silva and Remko Duursma
#'@references \link{http://maespa.github.io/}
#'@examples
#'
#'\dontrun{
#'
#'#' EXAMPLE 01: Plotting isolate trees
#'
#'#' shape cone crown 
#'open3d() 
#'PlotStand3D(crownshape = "cone", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'                        stemcolor = "chocolate4", resolution="high", shape=T) 
#'                        
#'#' elliptois crown shape 
#'open3d()
#'PlotStand3D(crownshape = "ellipsoid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'                        stemcolor = "chocolate4", resolution="high", shape=T) 
#'                        
#'#' halfellipsoid crown shape 
#'open3d()
#'PlotStand3D(crownshape = "halfellipsoid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'                        stemcolor = "chocolate4", resolution="high", shape=T) 
#'                        
#'#' paraboloid crown shape 
#'open3d()
#'PlotStand3D(crownshape = "paraboloid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'                        stemcolor = "chocolate4", resolution="high", shape=T) 
#'  
#'                                              
#'#'EXAMPLE 02: Plotting plantation forest stands
#' 
#'#' Set the lenght of the stand
#'xlenght<-30 # x lenght
#'ylenght<-20 # y lenght
#'
#'#' Set the space between trees
#'sx<-3 # x space lengh
#'sy<-2 # y space lenght
#'
#'#'#' Tree location grid
#'XYgrid <- expand.grid(x = seq(1,xlenght,sx), y = seq(1,ylenght,sy))
#'
#'#' Get the number of trees
#'Ntrees<-nrow(XYgrid)
#'
#'#'Plotting a Eucalyptus plantation stand using halfellipsoid for the tree crown shape
#'
#'#' Set stand trees parameters
#'meanHCB<-5 # mean tree crown base heigh
#'sdHCB<-0.1 # standard deviation tree crown base heigh
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # tree crown base heigh
#'CL<-HCB # tree crown heigh
#'CW<-HCB*0.6 # tree crown diameter
#'
#'open3d() # open a rgl window
#'
#'#'Plot stand
#'for( i in 1:Ntrees){
#'  PlotStand3D(crownshape = "halfellipsoid", CL = CL[i], CW = CW[i], 
#'              HCB = HCB[i], X = XYgrid[i,1], Y = XYgrid[i,2], dbh = 0.4, crowncolor = "forestgreen", 
#'                            stemcolor = "chocolate4", resolution="high", shape=T) 
#'                            }
#'                            
#'#' Plot parameters
#'axes3d(c("x-","x-", "y-","z"), col="gray") # axes
#'title3d(xlab = "Easting", ylab = "Northing", zlab = "Height", col="red") # title
#'planes3d(a=0,b=0,c=-1,d=0.0001,color="gray",alpha=1) # set a terrain plane
#'
#'
#'#'Plotting a Eucalyptus plantation stand using halfellipsoid for the tree crown shape
#'
#'#'Set stand trees parameters
#'meanHCB<-3 # mean tree crown base heigh
#'sdHCB<-0.1 # standard deviation tree crown base heigh
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # tree crown base heigh
#'CL<-HCB*2.0 # tree crown heigh
#'CW<-HCB*1.3 # tree crown diameter
#'
#'open3d() # open a rgl window
#'#' Plot stand
#'for( i in 1:Ntrees){
#'  PlotStand3D(crownshape = "cone", CL = CL[i], CW = CW[i], 
#'              HCB = HCB[i], X = XYgrid[i,1], Y = XYgrid[i,2], dbh = 0.4, crowncolor = "forestgreen", 
#'                            stemcolor = "chocolate4", resolution="high", shape=T) 
#'                            }
#'                            
#'#' Plot parameters
#'axes3d(c("x-","x-", "y-","z"), col="gray") # axes
#'title3d(xlab = "Easting", ylab = "Northing", zlab = "Height", col="red") # title
#'planes3d(a=0,b=0,c=-1,d=0.0001,color="gray",alpha=1) # set a terrain plane
#'
#'
#'#' EXAMPLE 03: Plotting natural mixed forest stands 
#'
#'#' Difers species of trees on the stand using diferents crown shapes
#'
#'#' Set the number of trees
#'Ntrees<-80 
#'
#'# Set the trees locations
#'xcoord<-sample(1:100,Ntrees) # x coord
#'coord<-sample(1:100,Ntrees) # x coord
#'
#'#'#' Set a location grid of trees 
#'XYgrid<-cbind(xcoord,ycoord)
#'
#'#' plot the location of the trees
#'plot(XYgrid, main="Tree location")
#'
#'meanHCB<-7 # mean tree crown base heigh
#'sdHCB<-3 # standard deviation tree crown base heigh
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # tree crown base heigh
#'crownshape<-sample(c("cone", "ellipsoid","halfellipsoid", "paraboloid"), Ntrees, replace=T) # tree crown shape 
#'CL<-HCB*1.3 # tree crown heigh
#'CW<-HCB # tree crown diameter
#'
#'open3d() # open a rgl window
#'#'Plot stand
#'
#'for( i in 1:Ntrees){
#'  PlotStand3D(crownshape = crownshape[i], CL = CL[i], CW = CW[i], 
#'              HCB = HCB[i], X = XYgrid[i,1], Y = XYgrid[i,2], dbh = 0.4, crowncolor = "forestgreen", 
#'                          stemcolor = "chocolate4", resolution="high", shape=T) 
#'                          }
#'                          
#'#' Plot parameters
#'axes3d(c("x-","x-", "y-","z"), col="gray") # axes
#'title3d(xlab = "Easting", ylab = "Northing", zlab = "Height", col="red") # title
#'planes3d(a=0,b=0,c=-1,d=0.0001,color="gray",alpha=1) # set a terrain plane
#'
#'
#'#' Difers trees height on the stand using diferents crown colors
#'
#'#'Set the number of trees
#'Ntrees<-80 
#'
#'#' Set the trees locations
#'xcoord<-sample(1:100,Ntrees) # x coord
#'ycoord<-sample(1:100,Ntrees) # x coord
#'
#'#'Set a location grid of trees 
#'XYgrid<-cbind(xcoord,ycoord)
#'
#'#'plot the location of the trees
#'plot(XYgrid, main="Tree location")
#'
#'meanHCB<-7 # mean tree crown base heigh
#'sdHCB<-3 # standard deviation tree crown base heigh
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # tree crown base heigh
#'crownshape<-sample(c("cone", "ellipsoid","halfellipsoid", "paraboloid"), Ntrees, replace=T) # tree crown shape 
#'CL<-HCB*1.3 # tree crown heigh
#'CW<-HCB # tree crown diameter
#'
#'#'Plot tree hiegh based on the HCB quantiles
#'HCBq<-quantile(HCB) # HCB quantiles
#'crowncolor<-NA*(1:Ntrees) # set a empty crowncolor vector
#'
#'#'classify trees by HCB quantile
#'for (i in 1:Ntrees){
#'  if (HCB[i] <= HCBq[2]) {crowncolor[i]<-"red"} # group 1
#'  if (HCB[i] > HCBq[2] & HCB[i] <= HCBq[3] ) {crowncolor[i]<-"blue"}  # group 2
#'  if (HCB[i] > HCBq[3] & HCB[i] <= HCBq[4] ) {crowncolor[i]<-"yellow"}  # group 3
#'  if (HCB[i] >= HCBq[4]) {crowncolor[i]<-"dark green"}  # group 4
#'  }
#'    
#'  open3d() # open a rgl window
#'#' Plot stand
#'for( i in 1:Ntrees){  
#'  PlotStand3D(crownshape = crownshape[i], CL = CL[i], CW = CW[i], 
#'    HCB = HCB[i], X = XYgrid[i,1], Y = XYgrid[i,2], dbh = 0.4, crowncolor = crowncolor[i], 
#'    stemcolor = "chocolate4", resolution="high", shape=T) 
#'    }
#'    
#'#' Plot parameters
#'axes3d(c("x-","x-", "y-","z"), col="gray") # axes
#'title3d(xlab = "Easting", ylab = "Northing", zlab = "Height", col="red") # title
#'planes3d(a=0,b=0,c=-1,d=0.0001,color="gray",alpha=1) # set a terrain plane
#'}
#' 
#'@export
LiDARtrees3D<-function (crownshape = c("cone", "ellipsoid",  
                                      "halfellipsoid", "paraboloid", "cylinder"), CL = 4, CW = 8, 
                       HCB = 10, X = 0, Y = 0, dbh = 0.3, crowncolor = "forestgreen", 
                       stemcolor = "chocolate4", resolution="high",shape=TRUE) 
{
  
  
  if (resolution=="low"){nz<-15;nalpha<-15}
  if (resolution=="median"){nz<-25;nalpha<-25}
  if (resolution=="high"){nz<-40;nalpha<-40}
  
  if ( shape==TRUE) {
    
  shape <- match.arg(crownshape)

  H <- HCB + CL
  dbase <- dbh * (H/(H - 1.3))
  if (!is.finite(dbase)) 
    dbase <- dbh
  
  
  m1 <- coord3dshape(shape, CW = CW, CL = CL, z0 = HCB, x0 = X, 
                                        y0 = Y, nz = nz, nalpha = nalpha)
  m2 <- coord3dshape("cone", CW = dbase, CL = H, z0 = 0, x0 = X, 
                      y0 = Y, nz = nz, nalpha = nalpha)
  
  interpol(m1, col = crowncolor, ...)
  interpol(m2, col = stemcolor, ...)
  
  } else {
    TreesModel(crownshape=crownshape, CW = CW, CL = CL, z0 = 0,HCB=HCB, x0 = X, 
                     y0 = Y, nz = nz, nalpha = nalpha, dbh = dbh,crowncolor = crowncolor, 
               stemcolor = stemcolor)
  }
  
  
}

coord3dshape <- function(crownshape=c("cone","ellipsoid","halfellipsoid","paraboloid","cylinder"),
                         nz=5, nalpha=5, CL=1, CW=1, x0=0, y0=0, z0=0
){

  crownshape <- match.arg(crownshape)
  
  z <- rep(seq(0,1,length=nz),each=nalpha)
  angs <- rep(seq(0,2*pi, length=nalpha),nz)

  if(crownshape == "cone")distfun <- (1-z)
  if(crownshape == "ellipsoid")distfun <- sqrt(1 - ((z-1/2)^2)/((1/2)^2))
  if(crownshape == "halfellipsoid")distfun <- sqrt(1 - z**2)
  if(crownshape == "paraboloid")distfun <- sqrt(1-z)
  if(crownshape == "cylinder")distfun <- 1
  
  
  r <- CW/2
  x <- x0 + r*distfun*cos(angs)
  y <- y0 + r*distfun*sin(angs)
  z <- z0 + z*CL
  
  keep <- !duplicated(cbind(x,y,z))
  x <- x[keep]
  y <- y[keep]
  z <- z[keep]
  return(matrix(cbind(x,y,z),ncol=3))
}

TreesModel<- function(crownshape=c("cone","ellipsoid","halfellipsoid","paraboloid","cylinder"),
                                   nz=5, nalpha=5, CL=5, CW=5, HCB=10, x0=0, y0=0, z0=0, dbh = 0.3, crowncolor = "forestgreen", 
                                      stemcolor = "chocolate4"
){
  
 crownshape <- match.arg(crownshape)
  
  z <- rep(seq(0,1,length=nz),each=nalpha)
  angs <- rep(seq(0,2*pi, length=nalpha),nz)

  if(crownshape == "cone")distfun <- (1-z)
  if(crownshape == "ellipsoid")distfun <- sqrt(1 - ((z-1/2)^2)/((1/2)^2))
  if(crownshape == "halfellipsoid")distfun <- sqrt(1 - z**2)
  if(crownshape == "paraboloid")distfun <- sqrt(1-z)
  if(crownshape == "cylinder")distfun <- 1
  H <- HCB + CL
  r <- CW/2
  x <- x0 + r*distfun*cos(angs)
  y <- y0 + r*distfun*sin(angs)
  z <- z0 + HCB + z*CL
  
  keep <- !duplicated(cbind(x,y,z))
  x <- x[keep]
  y <- y[keep]
  z <- z[keep]
  klj=matrix(cbind(x,y,z),ncol=3)
  
  mMatrix<-matrix(,ncol=3)[-1,]
  
  for ( i in 1:nrow(klj)){
    ln=i+nz
    
    if ( ln >= nrow(klj)) { ln2=nrow(klj) } else { ln2= ln}
    
    mMatrix<-rbind(mMatrix,rbind(klj[i,],klj[ln2,])) }
  
  
  kljzbase=subset(klj,klj[,3]==z[2])
  kljzbaseNew<-matrix(,ncol=3)[-1,]
  
  for ( i in 1:nrow(kljzbase)){
    kljzbaseNew<-rbind(kljzbaseNew,rbind(kljzbase[i,],c(x0,y0,HCB)))
    
  }
  
  newList<-rbind(kljzbaseNew,mMatrix,klj)
  plot3d(newList, type="l", col=crowncolor, add=T)
  m2 <- coord3dshape("cone", CW = dbh, CL = H, z0 = z0, x0 = x0, 
                     y0 = y0, nz = 50, nalpha = 50)
  interpol(m2, col = stemcolor)
  
}


interpol<- function(input,col) {
  surf.3d <- t(convhulln(input,options = "QJ")) 
  rgl.triangles(input[surf.3d,1],input[surf.3d,2],input[surf.3d,3],col=col,alpha = c(1.0),
                lit = TRUE,ambient = "black",specular = "white",emission = "black",shininess = 50.0,
                smooth = TRUE, texture = NULL,front = "fill",back ="fill",fog = F) 
}
