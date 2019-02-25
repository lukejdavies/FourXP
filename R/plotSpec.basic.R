#' Plots a basic spectrum
#'
#' @description Plots a basic spectrum with spec$flux and spec$wave paramters
#' 
#' @param spec spec object to plot
#' @param z FITS file row number  
#' @param degSmooth Hanning smoothing to apply (odd number)
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plotSpec.basic(spec=spec)
#' @export
plotSpec.basic<-function(spec=spec, cex.axis=1.4, cex.lab=1.4, degSmooth=7, xlim=NA, ylim=NA){
  
  options(warn=-1) 
  
  
  if (is.na(xlim)==T){xlim<-c(min(spec$wave, na.rm=T), max(spec$wave, na.rm=T))}
  if (is.na(ylim)==T){ylim<-c(min(spec$flux, na.rm=T), max(spec$flux, na.rm=T))}
  
  
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab=paste('Wavelength, ', spec$xunit,sep=''), ylab='Counts', grid=T, type='l', xlim=xlim, ylim=ylim, main=paste('ID=', spec$ID, ' - Hanning Smoothed, degree=',degSmooth,sep=''), cex.axis=1.4)

  options(warn=0)  
   
}