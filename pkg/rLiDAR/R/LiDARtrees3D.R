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
