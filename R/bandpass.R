#' Calculate spectral measurements within a specific filter bandbass
#'
#' @description Function convolves a given spectrum with a filter band pass and returns either
#' convolved spectrum or a total sum within the filter. Note that for the example to
#' work you must have unpacked the FourXP data files using:
#' 
#' @description > LibPaths<-.libPaths()
#' @description > system(paste('tar -xvf ', LibPaths, '/FourXP/data/calibrators.tar --directory ', LibPaths, '/FourXP/data/',sep='')) 
#' @description > system(paste('tar -xvf ', LibPaths, '/FourXP/data/idxFiles.tar --directory ', LibPaths, '/FourXP/data/' ,sep=''))
#' 
#' @param spec An R stucture containing spec$wave = vector of spectrum wavelengths
#' and spec$flux = vector of spectrum fluxes.     
#' @param filter A filter with which to convolve with. A 2d matrix with column
#' 1 as the filter wevelength (in the same units as spec$wave) and column
#' 2 the filter response. Could be prodcued by \code{getfilt}.
#' @param lum TRUE/FLASE, if TRUE \code{bandpass} will return total luminosity in the filter.
#' @author L. Davies, A. Robotham 
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plot(spec$wave, hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts')
#' filt<-getfilt('sdss_r')
#' lines(filt[,2],filt[,3]*max(spec$flux, na.rm=T), col='red')
#' lum<-bandpass(spec,cbind(filt[,2],filt[,3]), lum=T)
#' cat('Total luminosity in band = ',lum)
#' @export
bandpass=function(spec, filter, lum = T){
  wave<-spec$wave
  flux<-spec$flux
  tempfunc = approxfun(x = filter[, 1], y = abs(filter[, 2]))
  tempremap = tempfunc(wave)
  tempremap[is.na(tempremap)] = 0
  flux[is.na(flux)]= 0
  if (lum) {
    return = sum(tempremap * wave * flux)/sum(tempremap * wave)
  }
  else {
    return = tempremap * wave * flux/sum(tempremap * wave)
  }
}
