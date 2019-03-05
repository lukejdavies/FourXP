#' Load common spectral feature list
#'
#' @description Function to load common galaxy emission and absorption features. 
#' 
#' @return dataframe of line names, wavelengths in various units, frequency and emission/absorption flag  
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' #plot rest-frame spectrum
#' plot(spec$wave/(1+spec$z), hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts', ylim=c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*1.1))
#' lineList<-load.lines()
#' print(lineList)
#' lines(c(lineList$wave_ang[2], lineList$wave_ang[2]), c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*0.8), lty=2, col='blue')
#' text(lineList$wave_ang[2], max(spec$flux,na.rm=T),lineList$names[2], col='blue', cex=0.3)
#' lines(c(lineList$wave_ang[4], lineList$wave_ang[4]), c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*0.8), lty=2, col='blue')
#' text(lineList$wave_ang[4], max(spec$flux,na.rm=T),lineList$names[4], col='blue', cex=0.3)
#' lines(c(lineList$wave_ang[5], lineList$wave_ang[5]), c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*0.8), lty=2, col='blue')
#' text(lineList$wave_ang[5], max(spec$flux,na.rm=T),lineList$names[5], col='blue', cex=0.3)
#' lines(c(lineList$wave_ang[6], lineList$wave_ang[6]), c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*0.8), lty=2, col='blue')
#' text(lineList$wave_ang[6], max(spec$flux,na.rm=T),lineList$names[6], col='blue', cex=0.3)
#' lines(c(lineList$wave_ang[12], lineList$wave_ang[12]), c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*0.8), lty=2, col='darkgreen')
#' text(lineList$wave_ang[12], max(spec$flux,na.rm=T),lineList$names[12], col='darkgreen', cex=0.3)
#' lines(c(lineList$wave_ang[13], lineList$wave_ang[13]), c(-max(spec$flux,na.rm=T),max(spec$flux,na.rm=T)*0.8), lty=2, col='darkgreen')
#' text(lineList$wave_ang[13], max(spec$flux,na.rm=T),lineList$names[13], col='darkgreen', cex=0.3)
#' @export
load.lines=function(){
  c<-299792458

  names <- c('O VI','Ly-alpha', 'N V', 'O I','C II','Si IV','Si IV+O IV','C IV','He II', 'O III', 'Al III', 'C III','C II', 'Ne IV', 'Mg II', 'Ne V', 'Ne VI', 'O II', '', 'He I', 'S II','H-delta','H-gamma', 'O III', 'H-beta', 'O III', '','O III', 'O I', 'O I','N I','N II', 'H-alpha', 'N II', 'S II', 'S II', 'CaK&H', '', 'G-band', 'Mg', 'Na', 'Ca II', 'Ca II', 'Ca II') 
  wave_ang <- c(1033.82,1215.24,1240.81,1305.53,1335.31,1397.61, 1399.80, 1549.48, 1640.4, 1665.85, 1857.4,1908.734,2326.0,2439.5,2799.117,3346.79,3426.85,3727.092,3729.875,3889.0,4072.3,4102.89,4341.68,4364.436,4862.68,4932.603,4960.295,5008.240,6302.046,6365.536,6529.03,6549.86,6564.61,6585.27,6718.29,6732.67,3934.777,3969.588,4305.61,5176.7,5895.6,8500.36,8544.44,8664.52)
  stellar <- c(rep(FALSE,36), rep(TRUE,8))
  agnStrong <- c(T, T,T,F,F,F,T,T,F,F,F,T,T,F,T,F,F,F,F,F,T,T,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)
  strength <- c(1,9,3,0,0,0,1,8,0,0,0,7,0.5,0,8,0,0,5,0,0,0,2,3,0,4,0,2,3,0,0,0,3,8,3,3,3,1,1,1,1,1,1,1,1)+0.1
  freq_hz <- c/(wave_ang/(10.^10))
  wave_m <- wave_ang/(10.^10)
  wave_micron <- wave_ang/(10.^4)
  wave_nm <- wave_ang/(10.^1)
  lines<-data.frame(names, wave_ang,wave_m,wave_micron,wave_nm, freq_hz, stellar, agnStrong,strength)
  return(lines) 
}
