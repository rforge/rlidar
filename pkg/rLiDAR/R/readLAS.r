#'LiDAR - Reading LAS file
#'
#'@description Read LAS file format. The LAS file is a public file format for the interchange of LiDAR 3-dimensional point cloud data (American Society of Photogrammetry and Remote Sensing - ASPRS)
#'
#'@usage readLAS(LASfile, short=TRUE)
#'
#'@param LASfile A standard LAS data file accorgin to ASPRS
#'@param short Return only : x, y, z, intensity and return number
#'@return Return dataframe of the LAS data set
#'@author Michael Sumner and Carlos Alberto Silva
#'@examples
#'\dontrun{
#'# Importing LAS file:
#'myLAS<-data(LASfile) # or set a LAS  file (myLAS<-"LASfile.las")
#'
#'# Reading LAS file
#'readLAS(myLAS,short=TRUE)
#'
#'# Summary of the LAS file
#'summary(myLAS)
#'}
#'
#'@export
#'@importFrom bitops bitAnd bitShiftR
readLAS <- function(LASfile, short=TRUE) {

  skip <- 0
  nrows <- NULL
  
  hd <- publicHeaderDescription()
  pheader <- vector("list", nrow(hd))
  names(pheader) <- hd$Item
  con <- file(LASfile, open = "rb")
  isLASFbytes <- readBin(con, "raw", size = 1, n = 4, endian = "little")
  pheader[[hd$Item[1]]] <- readBin(isLASFbytes, "character", size = 4, endian = "little")
  if (! pheader[[hd$Item[1]]] == "LASF") {
    stop("Not a valid LAS file")
  }
  for (i in 2:nrow(hd)) {
    pheader[[hd$Item[i]]] <- readBin(con, what = hd$what[i], size = hd$Rsize[i], endian = "little", n = hd$n[i])
  }
  close(con)
  
  numberPointRecords <- pheader[["Number of point records"]]
  offsetToPointData <- pheader[["Offset to point data"]]
  pointDataRecordLength <-pheader[["Point Data Record Length"]]
  xyzScaleOffset <- cbind(unlist(pheader[c("X scale factor", "Y scale factor", "Z scale factor")]),
                          unlist(pheader[c("X offset", "Y offset", "Z offset")]))
  
  con <- file(LASfile, open = "rb")
  junk <- readBin(con, "raw", size = 1, n = offsetToPointData)
  
  if (skip > 0) {
    junk <- readBin(con, "raw", size = 1, n = pointDataRecordLength * skip)
    numberPointRecords <- numberPointRecords - skip
  }
  if (!is.null(nrows)) {
    if (numberPointRecords > nrows) numberPointRecords <- nrows
  }
  
  if (numberPointRecords < 1) stop("no records left to read")
  
  allbytes <- matrix(readBin(con, "raw", n = pointDataRecordLength * numberPointRecords, size = 1, endian = "little"),
                     ncol= pointDataRecordLength, nrow = numberPointRecords, byrow = TRUE)
  
  
  close(con)
  mm <- matrix(readBin(t(allbytes[,1:(3*4)]), "integer", size = 4, n = 3 * numberPointRecords, endian = "little"), ncol = 3, byrow = TRUE)
  
  mm[,1] <- mm[ ,1] * xyzScaleOffset[1,1] + xyzScaleOffset[1, 2]
  mm[,2] <- mm[ ,2] * xyzScaleOffset[2,1] + xyzScaleOffset[2, 2]
  mm[,3] <- mm[ ,3] * xyzScaleOffset[3,1] + xyzScaleOffset[3, 2]
  colnames(mm) <- c("X", "Y", "Z")
  
  Intensity <- readBin(t(allbytes[, 13:14]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    
  bytesList <- readBin(t(allbytes[,15]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")

  
  if (pheader[17][[1]]==00) {
    ReturnNumber <- bitAnd(7, bytesList)
    NumberOfReturns <- bitShiftR(bitAnd(56, bytesList), 3)
    ScanDirectionFlag <- bitShiftR(bitAnd(bytesList, 64), 6)
    EdgeofFlightLine <- bitShiftR(bitAnd(bytesList, 128), 7)
    Classification <- readBin(t(allbytes[, 16]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    ScanAngleRank <-readBin(t(allbytes[, 17]), "integer", size = 1, n = numberPointRecords, signed = TRUE, endian = "little")
    UserData <-readBin(t(allbytes[, 18]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    PointSourceID <-readBin(t(allbytes[, 19:20]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")

    if (short==TRUE) {cbind(mm, Intensity, ReturnNumber)} else {
      cbind(mm, Intensity, ReturnNumber,NumberOfReturns,ScanDirectionFlag,EdgeofFlightLine,Classification,ScanAngleRank,UserData,PointSourceID) }}
    
    
  if (pheader[17][[1]]==01) {
    ReturnNumber <- bitAnd(7, bytesList)
    NumberOfReturns <- bitShiftR(bitAnd(56, bytesList), 3)
    ScanDirectionFlag <- bitShiftR(bitAnd(bytesList, 64), 6)
    EdgeofFlightLine <- bitShiftR(bitAnd(bytesList, 128), 7)
    Classification <- readBin(t(allbytes[, 16]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    ScanAngleRank <-readBin(t(allbytes[, 17]), "integer", size = 1, n = numberPointRecords, signed = TRUE, endian = "little")
    UserData <-readBin(t(allbytes[, 18]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    PointSourceID <-readBin(t(allbytes[, 19:20]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    gpstime <- NULL
      
    if (ncol(allbytes) == 28) gpstime <- readBin(t(allbytes[ , 21:28]), "numeric", size = 8, n = numberPointRecords, endian = "little")
      
    if (short==TRUE) {cbind(mm, Intensity, ReturnNumber)} else {
        cbind(mm, Intensity, ReturnNumber,NumberOfReturns,ScanDirectionFlag,EdgeofFlightLine,Classification,ScanAngleRank,UserData,PointSourceID, gpstime) }}
      
        
  if (pheader[17][[1]]==02) {
    
    ReturnNumber <- bitAnd(7, bytesList)
    NumberOfReturns <- bitShiftR(bitAnd(56, bytesList), 3)
    ScanDirectionFlag <- bitShiftR(bitAnd(bytesList, 64), 6)
    EdgeofFlightLine <- bitShiftR(bitAnd(bytesList, 128), 7)
    Classification <- readBin(t(allbytes[, 16]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    
    ScanAngleRank <-readBin(t(allbytes[, 17]), "integer", size = 1, n = numberPointRecords, signed = TRUE, endian = "little")
    UserData <-readBin(t(allbytes[, 18]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    PointSourceID <-readBin(t(allbytes[, 19:20]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
 
    R <- readBin(t(allbytes[, 21:22]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    G <- readBin(t(allbytes[, 23:24]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    B <- readBin(t(allbytes[, 25:26]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
   
    
    if (short==TRUE) {cbind(mm, Intensity, ReturnNumber)} else {
      cbind(mm, Intensity, ReturnNumber,NumberOfReturns,ScanDirectionFlag,EdgeofFlightLine,Classification,ScanAngleRank,UserData,PointSourceID, R,G,B) }}
  
  if (pheader[17][[1]]==03) {
    
    ReturnNumber <- bitAnd(7, bytesList)
    NumberOfReturns <- bitShiftR(bitAnd(56, bytesList), 3)
    ScanDirectionFlag <- bitShiftR(bitAnd(bytesList, 64), 6)
    EdgeofFlightLine <- bitShiftR(bitAnd(bytesList, 128), 7)
    Classification <- readBin(t(allbytes[, 16]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    
    ScanAngleRank <-readBin(t(allbytes[, 17]), "integer", size = 1, n = numberPointRecords, signed = TRUE, endian = "little")
    UserData <-readBin(t(allbytes[, 18]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    PointSourceID <-readBin(t(allbytes[, 19:20]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    
    gpstime <- NULL
    gpstime <- readBin(t(allbytes[ , 21:28]), "numeric", size = 8, n = numberPointRecords, endian = "little")
    
    R <- readBin(t(allbytes[, 29:30]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    G <- readBin(t(allbytes[, 31:32]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    B <- readBin(t(allbytes[, 33:34]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    
    
    if (short==TRUE) {cbind(mm, Intensity, ReturnNumber)} else {
      cbind(mm, Intensity, ReturnNumber,NumberOfReturns,ScanDirectionFlag,EdgeofFlightLine,Classification,ScanAngleRank,UserData,PointSourceID, R,G,B) }}
  
  if (pheader[17][[1]]==04) {
    
    ReturnNumber <- bitAnd(7, bytesList)
    NumberOfReturns <- bitShiftR(bitAnd(56, bytesList), 3)
    ScanDirectionFlag <- bitShiftR(bitAnd(bytesList, 64), 6)
    EdgeofFlightLine <- bitShiftR(bitAnd(bytesList, 128), 7)
    Classification <- readBin(t(allbytes[, 16]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    
    ScanAngleRank <-readBin(t(allbytes[, 17]), "integer", size = 1, n = numberPointRecords, signed = TRUE, endian = "little")
    UserData <-readBin(t(allbytes[, 18]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    PointSourceID <-readBin(t(allbytes[, 19:20]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    
    gpstime <- NULL
    gpstime <- readBin(t(allbytes[ , 21:28]), "numeric", size = 8, n = numberPointRecords, endian = "little")
    
    WavePacket_Descriptor_Index <- readBin(t(allbytes[, 29]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    Byte_offset_to_waveform_data <- readBin(t(allbytes[, 30:37]), "integer", size = 8, n = numberPointRecords, signed = FALSE, endian = "little")
    Waveform_packet_size_in_bytes <- readBin(t(allbytes[, 38:41]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    Return_Point_Waveform_Location<- readBin(t(allbytes[, 42:45]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    X.t<- readBin(t(allbytes[, 46:49]), "intege", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    Y.t<- readBin(t(allbytes[, 50:53]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    Z.t<- readBin(t(allbytes[, 54:57]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    
    
    if (short==TRUE) {cbind(mm, Intensity, ReturnNumber)} else {
      cbind(mm, Intensity, ReturnNumber,NumberOfReturns,ScanDirectionFlag,EdgeofFlightLine,Classification,ScanAngleRank,UserData,PointSourceID,gpstime,
            WavePacket_Descriptor_Index,Byte_offset_to_waveform_data,Waveform_packet_size_in_bytes,
            Return_Point_Waveform_Location,X.t,Y.t,Z.t) }}
  
  
  if (pheader[17][[1]]==05) {
    
    ReturnNumber <- bitAnd(7, bytesList)
    NumberOfReturns <- bitShiftR(bitAnd(56, bytesList), 3)
    ScanDirectionFlag <- bitShiftR(bitAnd(bytesList, 64), 6)
    EdgeofFlightLine <- bitShiftR(bitAnd(bytesList, 128), 7)
    Classification <- readBin(t(allbytes[, 16]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    
    ScanAngleRank <-readBin(t(allbytes[, 17]), "integer", size = 1, n = numberPointRecords, signed = TRUE, endian = "little")
    UserData <-readBin(t(allbytes[, 18]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    PointSourceID <-readBin(t(allbytes[, 19:20]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    
    gpstime <- NULL
    gpstime <- readBin(t(allbytes[ , 21:28]), "numeric", size = 8, n = numberPointRecords, endian = "little")
    
    R <- readBin(t(allbytes[, 29:30]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    G <- readBin(t(allbytes[, 31:32]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    B <- readBin(t(allbytes[, 33:34]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    
    WavePacket_Descriptor_Index <- readBin(t(allbytes[, 35]), "integer", size = 1, n = numberPointRecords, signed = FALSE, endian = "little")
    Byte_offset_to_waveform_data <- readBin(t(allbytes[, 36:43]), "integer", size = 8, n = numberPointRecords, signed = FALSE, endian = "little")
    Waveform_packet_size_in_bytes <- readBin(t(allbytes[, 44:47]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    Return_Point_Waveform_Location<- readBin(t(allbytes[, 48:51]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    X.t<- readBin(t(allbytes[, 52:55]), "intege", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    Y.t<- readBin(t(allbytes[, 56:59]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    Z.t<- readBin(t(allbytes[, 60:63]), "integer", size = 4, n = numberPointRecords, signed = FALSE, endian = "little")
    
    if (short==TRUE) {cbind(mm, Intensity, ReturnNumber)} else {
      cbind(mm, Intensity, ReturnNumber,NumberOfReturns,ScanDirectionFlag,EdgeofFlightLine,Classification,ScanAngleRank,UserData,PointSourceID,gpstime,
            R, G, B, WavePacket_Descriptor_Index,Byte_offset_to_waveform_data,Waveform_packet_size_in_bytes,
            Return_Point_Waveform_Location,X.t,Y.t,Z.t) }}
 
}


publicHeaderDescription <- function() {
  hd <- structure(list(Item = c("File Signature (\"LASF\")",
                                "(1.1) File Source ID", "(1.1) Global Encoding",
                                "(1.1) Project ID - GUID data 1", "(1.1) Project ID - GUID data 2",
                                "(1.1) Project ID - GUID data 3", "(1.1) Project ID - GUID data 4",
                                "Version Major", "Version Minor", "(1.1) System Identifier",
                                "Generating Software", "(1.1) File Creation Day of Year",
                                "(1.1) File Creation Year", "Header Size", "Offset to point data",
                                "Number of variable length records",
                                "Point Data Format ID (0-99 for spec)", "Point Data Record Length",
                                "Number of point records", "Number of points by return",
                                "X scale factor", "Y scale factor", "Z scale factor", "X offset",
                                "Y offset", "Z offset", "Max X", "Min X", "Max Y", "Min Y", "Max Z",
                                "Min Z"), Format = c("char[4]", "unsigned short", "unsigned short",
                                                     "unsigned long", "unsigned short", "unsigned short",
                                                     "unsigned char[8]", "unsigned char", "unsigned char", "char[32]",
                                                     "char[32]", "unsigned short", "unsigned short", "unsigned short",
                                                     "unsigned long", "unsigned long", "unsigned char", "unsigned short",
                                                     "unsigned long", "unsigned long[5]", "double", "double", "double",
                                                     "double", "double", "double", "double", "double", "double", "double",
                                                     "double", "double"), Size = c("4 bytes", "2 bytes", "2 bytes",
                                                                                   "4 bytes", "2 byte", "2 byte", "8 bytes", "1 byte", "1 byte",
                                                                                   "32 bytes", "32 bytes", "2 bytes", "2 bytes", "2 bytes", "4 bytes",
                                                                                   "4 bytes", "1 byte", "2 bytes", "4 bytes", "20 bytes", "8 bytes",
                                                                                   "8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes",
                                                                                   "8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes"), Required =
                         c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*",
                           "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*",
                           "*", "*", "*", "*", "*")), .Names = c("Item", "Format", "Size",
                                                                 "Required"), row.names = 2:33, class = "data.frame")
  hd$what <- ""
  hd$what[grep("unsigned", hd$Format)] <- "integer"
  hd$what[grep("char", hd$Format)] <- "raw"
  hd$what[grep("short", hd$Format)] <- "integer"
  hd$what[grep("long", hd$Format)] <- "integer"
  hd$what[grep("double", hd$Format)] <- "numeric"
  hd$signed <- TRUE
  hd$signed[grep("unsigned", hd$Format)] <- FALSE
  hd$n <- as.numeric(gsub("[[:alpha:][:punct:]]", "", hd$Format))
  hd$n[hd$what == "character"] <- 1
  hd$n[is.na(hd$n)] <- 1
  hd$Hsize <- as.numeric(gsub("[[:alpha:]]", "", hd$Size))
  hd$Rsize <- hd$Hsize / hd$n
  hd$Rsize[hd$what == "raw"] <- 1
  hd$n[hd$what == "raw"] <- hd$Hsize[hd$what == "raw"]
  hd
}
