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
  spec$z<-NA
  spec$expMin<-NA
  spec$prob<-NA
  
  spec2<-spec
  spec2$wave<-seq(min(wave), max(wave),0.5)
  spec2$flux<-approx(wave, flux, spec2$wave)$y
  
  spec3<-scFlux(spec2, spec$mag, "sdss_r")
  sc<-median(spec3$flux)/median(spec$flux)
  spec3$sky<-spec3$sky/sc
  spec3$var<-spec3$var/sc
  spec3$error<-spec3$error/sc
  
  return(spec3)

}

