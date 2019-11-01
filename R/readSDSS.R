#' Read SDSS spectra
#'
#' @description Reads in an SDSS sepctrum and converts to a sensibele format for 4XP
#' @param file SDSS spectrum filename to read   
#' @param id id to assing to spectum

#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 

#' @export
readSDSS<-function(file, id='test'){
  specA<-readFITS(file)
  data<-readFITS(file, hdu=2)
  spec<-list()
  flux<-specA$col[1][[1]]*10^-17
  wave<-10^specA$col[2][[1]]
  spec$id<-id
  spec$wave<-wave
  spec$flux<-flux
  spec$xunit<-'ang'
  spec$yunit<-'ang'
  spec$sky<-specA$col[7][[1]]*10^-17
  spec$var<-(1/(specA$col[3][[1]]))*10^-17
  spec$error<-sqrt(1/specA$col[3][[1]])*10^-17
  spec$mag<-magABspec(spec,filter = "sdss_r")
  spec$band<-"sdss_r"
  spec$col<-NA
  spec$mass<-NA
  spec$sfr<-NA
  spec$agn<-NA
  spec$expMin<-NA
  spec$prob<-NA
  
  spec2<-spec
  spec2$wave<-seq(min(wave), max(wave),0.5)
  spec2$flux<-approx(wave, flux, spec2$wave)$y
  spec2$sky<-approx(wave, spec$sky, spec2$wave)$y
  spec2$var<-approx(wave, spec$var, spec2$wave)$y
  spec2$error<-approx(wave, spec$error, spec2$wave)$y
  
  spec3<-scFlux(spec2, spec$mag, "sdss_r")
  sc<-median(spec3$flux)/median(spec$flux)
  spec3$sky<-spec3$sky/sc
  spec3$var<-spec3$var/sc
  spec3$error<-spec3$error/sc
  
  
  specDat<-array(NA, dim=c(1,length(data$colNames)))
  for (i in 1:length(data$colNames)){
    specDat[1,i]<-data$col[[i]][1]
  }
  colnames(specDat)<-data$colNames
  spec3$specDat<-specDat
  
  spec3$z<-as.numeric(spec3$specDat[,'Z'])
  spec3$RA<-as.numeric(spec3$specDat[,'PLUG_RA'])
  spec3$DEC<-as.numeric(spec3$specDat[,'PLUG_DEC'])
  spec3$UTMJD<-as.double(2400000.0) + as.double(spec3$specDat[,'MJD']) / (24.*3600.)
  spec3$longitude <- -105.820417
  spec3$latitude <- 32.780361
  spec3$altitude<-2788
  
  
  
  return(spec3)

}

