#' Read in spectra from FITS files 
#'
#' @description Reads in spectra from FITS files in a sensible-ish way and make list useable by
#'  other functions in FourXP. 
#' 
#' @param file file name of FITS file to read
#' @param row row number to read in
#' @param logw TRUE/FALSE is file in log wavelength
#' @param sc value to scale inout fluxes
#' @param xunit xunit value of spectrum (ang, hz, micron,nm,m)
#' @param yunit xunit value of spectrum (ang=ergs/sec/cm^2/ang, hz=ergs/sec/cm^2/hz, Jy)
#' @param z input redshfit for spectrum (tries to read this from header)
#' @param RA input RA for spectrum (tries to read this from header)
#' @param DEC input DEC for spectrum (tries to read this from header)
#' @param errorRow row in spectrum of error (if available)
#' @param SDSS_spec TRUE/FALSE is this an SDSS spectrum (uses specific format) 
#' @return List containing wave,flux,xunit,yunit, z, RA, DEC, error
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' name<-paste(specDir,'Spectra/AGN_v1.0_z050_mr165_type1.fits',sep='')
#' spec<-get.spec(name, xunit='ang', yunit='ang', z=0.5)
#' plotSpec(spec, Air=F)
#' @export
get.spec=function(file, row=1, logw=F, sc=1, xunit='ang', yunit='ang', z=NA, RA=NA, DEC=NA, errorRow=NA, SDSS_spec=F) { 
    specf<-readFITS(file)

    if (SDSS_spec==F){


        if (length(as.numeric(specf$hdr[which(specf$hdr=='Z')+1])) >0 & is.finite(z)==F) {z<-as.numeric(specf$hdr[which(specf$hdr=='Z')+1])}
        if (length(as.numeric(specf$hdr[which(specf$hdr=='RA')+1])) >0 & is.finite(RA)==F) {RA<-as.numeric(specf$hdr[which(specf$hdr=='RA')+1])}
        if (length(as.numeric(specf$hdr[which(specf$hdr=='DEC')+1])) >0 & is.finite(DEC)==F) {DEC<-as.numeric(specf$hdr[which(specf$hdr=='DEC')+1])}
        if (length(dim(specf$imDat))>1) {flux <- specf$imDat[,row]*sc}
        if (length(dim(specf$imDat))==1) {flux <- specf$imDat*sc}
        
        if (is.finite(errorRow)==T) {
            if(dim(specf$imDat)[2]>=errorRow){error <- specf$imDat[,errorRow]*sc} 
            else {
                error<-flux
                error[]<-NA
                cat(' ** WARNING ERROR ROW NOT AVAILBLE, NO ERROR OUTPUT **', '\n')
            }
        }else{
            error<-flux
            error[]<-NA
        }
        
        CRVAL<-as.numeric(specf$hdr[which(specf$hdr=='CRVAL1')+1])
        CRPIX<-as.numeric(specf$hdr[which(specf$hdr=='CRPIX1')+1]) 
        CDELT<-as.numeric(specf$hdr[which(specf$hdr=='CDELT1')+1])
        if (length(CDELT)==0) {CDELT<-as.numeric(specf$hdr[which(specf$hdr=='CD1_1')+1])}
        
        wave <- ((CRVAL+(c(0:(specf$axDat$len[1]-1))*CDELT))-(CDELT*(CRPIX-1)))
        if (logw==T) {
            wave<-wave-CDELT
            wave<-10.^wave
        }
    }

    if (SDSS_spec==T){

        wave<-10.^specf$col[[2]]
        flux<-specf$col[[1]]
        error<-specf$col[[3]]

        }

    
  spec <- list(wave,flux,xunit,yunit, z, RA, DEC, error)
  names(spec) <- c('wave', 'flux', 'xunit', 'yunit', 'z','RA','DEC', 'error')
  return=spec      
}

