#' Plots a basic spectrum
#'
#' @description Plots a basic spectrum with spec$flux and spec$wave parameters
#' 
#' @param spec spec object to plot
#' @param z FITS file row number  
#' @param degSmooth Hanning smoothing to apply (odd number)
#' @param xlim x range of plot (if not supplied will use sensible values)
#' @param ylim y range of plot (if not supplied will use sensible values)
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plotSpec.basic(spec=spec)
#' @export
plotSpec.basic<-function(spec=spec, cex.axis=1.4, cex.lab=1.4, degSmooth=7, xlim=NA, ylim=NA){
  
  options(warn=-1) 
  
  
  if (is.na(xlim)==T){xlim<-c(min(spec$wave, na.rm=T), max(spec$wave, na.rm=T))}
  if (is.na(ylim)==T){ylim<-c(min(spec$flux, na.rm=T), max(hanning.smooth(spec$flux, degree=degSmooth), na.rm=T))}
  
  
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab=paste('Wavelength, ', spec$xunit,sep=''), ylab='Counts', grid=T, type='l', xlim=xlim, ylim=ylim, main=paste('ID=', spec$id, ' - Hanning Smoothed, degree=',degSmooth,sep=''), cex.axis=1.4)

  options(warn=0)  
   
}