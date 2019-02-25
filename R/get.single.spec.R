#' Extracts single spectrum from FITS frame row 
#'
#' @description Extracts a spectrum from a given row in a FITS file
#' 
#' @param file The file location 
#' @param SpecNum FITS file row number  
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @param verbose 0=quiet, 1=tell me some things that are going on,  2=tell me all things that are going on
#' @examples 

#' @export
get.single.spec<-function(file=file, hdu=1, SpecNum=1, tryMeta=FALSE, xunit='Angstrom', yunit='Counts', ID=NA, verbose=2){
  
  options(warn=-1) 
  
  im <- readFITS(file=file,hdu=hdu)
  CRVAL1 <- as.numeric(im$hdr[which(im$hdr=="CRVAL1")+1])
  CRPIX1 <- as.numeric(im$hdr[which(im$hdr=="CRPIX1")+1])
  CDELT1 <- as.numeric(im$hdr[which(im$hdr=="CDELT1")+1])
  wave <- CRVAL1+((c(1:dim(im$imDat)[1])-CRPIX1)*CDELT1)
  flux<-im$imDat[,SpecNum]

  RA<-NA
  DEC<-NA
  z<-NA
  
  if (tryMeta==T){
    
    if (length(im$hdr[which(im$hdr=="OBJECT")+1])>0) {ID<-im$hdr[which(im$hdr=="OBJECT")+1][1]}
    if (length(im$hdr[which(im$hdr=="RA")+1])>0) {RA<-as.numeric(im$hdr[which(im$hdr=="RA")+1])}
    if (length(im$hdr[which(im$hdr=="DEC")+1])>0) {DEC<-as.numeric(im$hdr[which(im$hdr=="DEC")+1])}
    if (length(im$hdr[which(im$hdr=="Z")+1])>0) {z<-as.numeric(im$hdr[which(im$hdr=="Z")+1])}
    if (length(im$hdr[which(im$hdr=="CUNIT1")+1])>0) {xunit<-(im$hdr[which(im$hdr=="CUNIT1")+1])}
  }
  
  if (is.na(ID)==T){ID<-SpecNum}
  
  spec<-list(ID=ID, wave=wave,flux=flux, RA=RA, DEC=DEC, z=z, xunit=xunit, yunit=yunit)

  options(warn=0)
  
  return(spec)
  }