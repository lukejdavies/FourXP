#' Extracts spectra from AAT+AAOmega FITS frame 
#'
#' @description Takes a reduced AAOOmegam frame and extracts spectra numbers given by
#' a vector, SpecNums. Produces a list with one object per spectrum containong all meta information 
#' wavelength and flux.
#' 
#' @param file The file location of the AAOmega frame
#' @param SpecNums SpectraNumbers to Extract   
#' @param verbose 0=quiet, 1=tell me some things that are going on,  2=tell me all things that are going on
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' allSpec<-load.spec(file="2018_1_15_config_6_reduced.fits", SpecNums=c(23,45,89), verbose=2)
#' names(allSpec)
#' plotSpec(spec=allSpec$spec23)
#' @export
load.spec<-function(file=file, SpecNums=1, verbose=2){

  if (verbose>0){cat('- Reading FITS file....', '\n')}
  
tab <- readFITS(file=file,hdu=3)
im <- readFITS(file=file,hdu=1)
sn <- readFITS(file=file,hdu=2)
sky <- readFITS(file=file,hdu=8)
CRVAL1 <- as.numeric(im$hdr[which(im$hdr=="CRVAL1")+1])
CRPIX1 <- as.numeric(im$hdr[which(im$hdr=="CRPIX1")+1])
CDELT1 <- as.numeric(im$hdr[which(im$hdr=="CDELT1")+1])
wave <- CRVAL1+((c(1:dim(im$imDat)[1])-CRPIX1)*CDELT1)
UTMJD<-as.double(tab$hdr[which(tab$hdr=='CONFMJD')+1])
CONFIG<-im$hdr[which(im$hdr=="CFG_FILE")+1]
EXP<-as.numeric(im$hdr[which(im$hdr=="EXPOSED")+1])
allSpec<-list()

if (verbose>0){cat('- Looping over spectra....', '\n')}

for (j in SpecNums) {

  ID<-strsplit(tab$col[[which(tab$colNames=='NAME')]][j],' ')[[1]][1]

  if (verbose>1){cat('      - Extracting spectrum ', j, ' of ',dim(im$imDat)[2],': ',ID, '\n')}
  
  RA<-tab$col[which(tab$colNames=='RA')][[1]][j]*(180/pi)
  DEC<-tab$col[which(tab$colNames=='DEC')][[1]][j]*(180/pi)
  X<-tab$col[which(tab$colNames=='X')][[1]][j]
  Y<-tab$col[which(tab$colNames=='Y')][[1]][j]
  XERR<-tab$col[which(tab$colNames=='XERR')][[1]][j]
  YERR<-tab$col[which(tab$colNames=='YERR')][[1]][j]
  THETA<-tab$col[which(tab$colNames=='THETA')][[1]][j]
  TYPE<-tab$col[which(tab$colNames=='TYPE')][[1]][j]
  
  if (TYPE=='P' & substr(ID,1,1)!='S'){TYPE<-'target'}
  if (TYPE=='P' & substr(ID,1,1)=='S'){TYPE<-'std'}
  if (TYPE=='F'){TYPE<-'guide'}
  if (TYPE=='S'){TYPE<-'sky'}
  
  
  PIVOT<-tab$col[which(tab$colNames=='PIVOT')][[1]][j]
  MAG<-tab$col[which(tab$colNames=='MAGNITUDE')][[1]][j]
  FIBRE<-j
  
  fluxSpec<-im$imDat[,j]
  snSpec <- sn$imDat[,j]
  skySpec <- sky$imDat
  
  
  specTmp<-list(wave=wave,flux=fluxSpec,sn=snSpec,sky=skySpec, ID=ID, origfile=file, RA=RA, DEC=DEC, X=X, Y=Y,XERR=XERR,YERR=YERR,THETA=THETA, TYPE=TYPE,PIVOT=PIVOT,MAG=MAG, FIBRE=FIBRE, UTMJD=UTMJD,DATE=date, CONFIG=CONFIG, xunit='ang', yunit='counts', z=NA, EXP=EXP)
  assign(paste('spec',j,sep=''), specTmp)

  allSpec[[paste('spec',j,sep='')]] <- specTmp
  
}

return(allSpec)
}