#' Measure Spectral AB magnitude in given filter
#'
#' @description Convolves an input spectrum with a user-defined filter and calculate 
#' the spectral AB magnitude. Available filters are given in getfilt: GALEX_FUV, GALEX_NUV, 
#' sdss_u,sdss_g,sdss_r,sdss_i,sdss_z, des_u,des_g, des_r,des_i,des_z,des_Y,VISTA_Z,VISTA_Y,
#' VISTA_J,VISTA_H,VISTA_K, WISE_W1, WISE_W2,WISE_W3,WISE_W4,Herschel_100,Herschel_160,
#' Herschel_250,Herschel_350,Herschel_500. Note that for this to work you must have unpacked 
#' the TAZ data files using:
#' 
#' @description > LibPaths<-.libPaths()
#' @description > system(paste('tar -xvf ', LibPaths, '/FourXP/data/calibrators.tar --directory ', LibPaths, '/FourXP/data/',sep='')) 

#'  
#' @param spec An R stucture containing spec$wave = vector of spectrum wavelengths
#' and spec$flux = vector of spectrum fluxes, spec$xunit = Units in the x-direction 
#' (Allowed units are 'ang', 'hz', 'm', 'micron','nm'), spec$yunit = Units in the y-direction 
#' (Allowed units are 'ang'=ergs/sec/cm^2/Ang, 'hz'=ergs/sec/cm^2/Hz, 'Jy'=Jansky.)
#' @param filter string of the filter to load.
#' @return AB magnituide of spectrum in given filter  
#' @author A. Robotham
#' @author L. Davies
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' 
#' #scale to ergs/sec/cm^2/ang
#' spec$flux<-spec$flux*spec$fluxSc
#' 
#' plot(spec$wave, hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts')
#' 
#' filt<-getfilt('sdss_r')
#' lines(filt[,2],filt[,3]*max(spec$flux, na.rm=T), col='red')
#' 
#' print(magABspec(spec, filter='sdss_r'))

#' @export
magABspec=function(spec, filter='sdss_r'){
  if (spec$xunit[1]=='ang') {wavefac=1e-10}
  if (spec$xunit[1]=='mircon') {wavefac=1e-6}
  if (spec$xunit[1]=='m') {wavefac=1}
  if (spec$xunit[1]=='nm') {wavefac=1e-9} 
  c<-299792458

  wave<-spec$wave
  flux<-spec$flux

  mag<-c(1:length(filter))
  for (i in 1:length(filter)){
    
    filter2=getfilt(filter[i])[,2:3]

    if (spec$yunit[1]=='ang') {fluxnu=(wavefac*flux*wave^2)/c}
    if (spec$yunit[1]=='hz') {fluxnu=flux}
    if (spec$yunit[1]=='Jy') {fluxnu=flux*1e-23}
    
    spec$flux<-fluxnu
    
    totlumnu = bandpass(spec, filter = filter2, lum = T)
    
    mag[i]= -2.5 * log10(totlumnu) - 48.6
  }
  return=mag

}

