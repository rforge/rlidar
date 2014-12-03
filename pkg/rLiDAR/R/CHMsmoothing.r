#'@export
CHMsmoothing<-function(chm, filter="mean", ws=5, sigma=0.6) {

  if (class(chm)[1]!='RasterLayer') {
      chmInput<-raster(chm)
      } else {chmInput<-chm
  }
  
  if (filter == "mean") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=mean)
  }
  if (filter == "median") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=median)
  }
  if (filter == "maximum") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=max)
  }
  if (filter == "minimum") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=min)
  }
  
  if (filter =="gaussian") {
    
    fgauss <- function(sigma, n=ws) {
      m <- matrix(nc=n, nr=n)
      col <- rep(1:n, n)
      row <- rep(1:n, each=n)
      x <- col - ceiling(n/2)
      y <- row - ceiling(n/2)
      m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
      m / sum(m)
    }
    gf=fgauss(sigma)
    chmR <- focal(chmInput, w=gf)
  }
  
  return(chmR)
}
