\name{chm}
\alias{chm}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(chm)}
\format{
  The format is:
Formal class 'RasterLayer' [package "raster"] with 12 slots
  ..@ file    :Formal class '.RasterFile' [package "raster"] with 13 slots
  .. .. ..@ name        : chr ""
  .. .. ..@ datanotation: chr "FLT4S"
  .. .. ..@ byteorder   : chr "little"
  .. .. ..@ nodatavalue : num -Inf
  .. .. ..@ NAchanged   : logi FALSE
  .. .. ..@ nbands      : int 1
  .. .. ..@ bandorder   : chr "BIL"
  .. .. ..@ offset      : int 0
  .. .. ..@ toptobottom : logi TRUE
  .. .. ..@ blockrows   : int 0
  .. .. ..@ blockcols   : int 0
  .. .. ..@ driver      : chr ""
  .. .. ..@ open        : logi FALSE
  ..@ data    :Formal class '.SingleLayerData' [package "raster"] with 13 slots
  .. .. ..@ values    : num [1:25281] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ offset    : num 0
  .. .. ..@ gain      : num 1
  .. .. ..@ inmemory  : logi TRUE
  .. .. ..@ fromdisk  : logi FALSE
  .. .. ..@ isfactor  : logi FALSE
  .. .. ..@ attributes: list()
  .. .. ..@ haveminmax: logi TRUE
  .. .. ..@ min       : num 0
  .. .. ..@ max       : num 28.3
  .. .. ..@ band      : int 1
  .. .. ..@ unit      : chr ""
  .. .. ..@ names     : chr "band1"
  ..@ legend  :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. ..@ type      : chr(0) 
  .. .. ..@ values    : logi(0) 
  .. .. ..@ color     : logi(0) 
  .. .. ..@ names     : logi(0) 
  .. .. ..@ colortable: logi(0) 
  ..@ title   : chr(0) 
  ..@ extent  :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. ..@ xmin: num 747598
  .. .. ..@ xmax: num 747678
  .. .. ..@ ymin: num 3460081
  .. .. ..@ ymax: num 3460161
  ..@ rotated : logi FALSE
  ..@ rotation:Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. ..@ geotrans: num(0) 
  .. .. ..@ transfun:function ()  
  ..@ ncols   : int 159
  ..@ nrows   : int 159
  ..@ crs     :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. ..@ projargs: chr NA
  ..@ history : list()
  ..@ z       : list()
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(chm)
## maybe str(chm) ; plot(chm) ...
}
\keyword{datasets}
