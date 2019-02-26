#' Plots a spectrum object generated as part of TAZ
#'
#' @description Function takes a 'spec' R list with sepc$wave, spec$flux, ...
#' and plots the full 1D spectrum and zoom ins of key emission and absorption 
#' lines. Lines will be scaled to the correct redshfit in spec$z is provided. If 
#' no spec$z is provided, line positions and zooms are plotted in the rest-frame. 
#' 
#' @param spec An R stucture containing spec$wave = vector of spectrum wavelengths
#' and spec$flux = vector of spectrum fluxes.  
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @param leg.cex Size of font in legend   
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plotSpec(spec)
#' @export
plotSpec<-function(spec=spec, leg.cex=1.4, degSmooth=7, xlim=NA, ylim=NA){
  
  
  

  
  if (is.na(xlim[1])==T){xlim<-c(min(spec$wave, na.rm=T), max(spec$wave, na.rm=T))}
  if (is.na(ylim[1])==T){ylim<-c(min(spec$flux, na.rm=T), max(spec$flux, na.rm=T))}
  
options(warn=-1)  
  

  par(mfrow = c(3, 3))
  par(mar=c(3.1,3.1,1.1,1.1))
  
  layout(matrix(c(1,1,1, 2,3,4, 5,6,7), 3, 3, byrow = TRUE))
  
  
  textZ<-spec$z
  
  if (is.null(spec$MAG)==TRUE){spec$MAG<-NA}
  if (is.null(spec$ID)==TRUE){spec$ID<-'No ID Provided'}  
  if (is.null(spec$EXP)==TRUE){spec$EXP<-NA}
  if (is.null(spec$prob)==TRUE){spec$prob<-NA} 
  if (is.null(spec$prob)==TRUE){spec$prob<-NA} 
  
  spec$EXP<-as.numeric(spec$EXP)
  
  if (is.null(spec$z)==TRUE){
    spec$z<-0
    textZ<-NA
  }
  
  
  
  lines<-load.lines()
  line_x <- as.numeric(lines$wave_ang)*(1+spec$z)
  
  plotZ<-spec$z
  if (is.finite(plotZ)==FALSE){
    plotZ<-0
    textZ<-NA
  }
  
  
  peak_F<-max(spec$flux, na.rm=T)
  
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=xlim, ylim=ylim, main=paste(spec$ID, ' - Hanning Smoothed, degree=',degSmooth,sep=''), cex.axis=1.4)
  
  
  
  plotLines(z=plotZ, xunit='ang', labPos=0.8*max(spec$flux,na.rm=T), lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-40)
  
  if (is.null(spec$tmpWave)==TRUE){spec$tmpWave<-spec$wave}
  if (is.null(spec$tmpFlux)==TRUE){spec$tmpFlux<-rep(NA, length(spec$wave))}
  if (is.null(spec$sky)==TRUE){spec$sky<-rep(NA, length(spec$wave))}
  if (is.null(spec$sn)==TRUE){spec$sn<-rep(NA, length(spec$wave))}
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(spec$flux,na.rm=T)-min(spec$flux,na.rm=T))*0.4)-((ylim[2]-ylim[1])*0.3)
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(spec$flux,na.rm=T)-min(spec$flux,na.rm=T))*0.4)-((ylim[2]-ylim[1])*0.45)
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(spec$flux,na.rm=T)-min(spec$flux,na.rm=T))*0.4)-((ylim[2]-ylim[1])*0.5)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  
  
  
  legend('topleft', legend=c(paste('ID=',spec$ID,sep=''), paste('z=',format(textZ, nsmall=2,digits=2),sep=''), paste('Y-mag=',spec$MAG,sep=''), paste('Prob=',format(spec$prob, nsmall=2,digits=2),sep=''), paste('TEXP (min)=',spec$EXP/60,sep=''), 'Sky', 'Variance', 'BestFit Template'), text.col=c('black','black','black','black','black','darkgreen', 'gold', 'indianred2'), bg='white', cex=leg.cex)
  
  
  
  textZ<-spec$z
  
  if (is.null(spec$MAG)==TRUE){spec$MAG<-NA}
  if (is.null(spec$ID)==TRUE){spec$ID<-'No ID Provided'}  
  if (is.null(spec$EXP)==TRUE){spec$EXP<-NA}
  if (is.null(spec$prob)==TRUE){spec$prob<-NA} 
  
  if (is.null(spec$z)==TRUE){
    spec$z<-0
    textZ<-NA
  }
  
  
  
  lines<-load.lines()
  line_x <- as.numeric(lines$wave_ang)*(1+spec$z)
  
  plotZ<-spec$z
  if (is.finite(plotZ)==FALSE){
    plotZ<-0
    textZ<-NA
  }
  
  
  peak_F<-max(spec$flux, na.rm=T)
  
  
  
  range<-250
  wavePoint<-3727
  peak<-1.2*max(spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T)
  if (is.finite(peak)==F){peak=peak_F}
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=(1+plotZ)*wavePoint+c(-range,range), ylim=c(-peak,peak), main='OII', cex.axis=1.4)
  #plotLines(z=plotZ, xunit='ang', labPos=0.8*peak, lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-20)
  
  
  lines((1+plotZ)*c(3728, 3728), c(-100,100), lty=2, col='blue')
  text(((1+plotZ)*c(3728))-20, 0.8*peak, 'OII', srt=270, col='blue', cex=1.2)
  
  yr<-c(-peak,peak)
  hereFlux<-spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)]
  
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sky<-sky+abs(median(sky[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.45)
  
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sn<-sn+abs(median(sn[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.5)
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.4)
  tmp<-tmp+abs(median(tmp[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.3)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  range<-500
  wavePoint<-4950
  peak<-1.2*max(spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T)
  if (is.finite(peak)==F){peak=peak_F}
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=(1+plotZ)*wavePoint+c(-range,range), ylim=c(-peak,peak), main='OIII/H-beta', cex.axis=1.4)
  plotLines(z=plotZ, xunit='ang', labPos=0.8*peak, lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-20)
  yr<-c(-peak,peak)
  hereFlux<-spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)]
  
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sky<-sky+abs(median(sky[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.45)
  
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sn<-sn+abs(median(sn[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.5)
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.4)
  tmp<-tmp+abs(median(tmp[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.3)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  range<-300
  wavePoint<-6650
  peak<-1.2*max(spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T)
  if (is.finite(peak)==F){peak=peak_F}
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=(1+plotZ)*wavePoint+c(-range,range), ylim=c(-peak,peak), main='H-alpha, NII, SII', cex.axis=1.4)
  plotLines(z=plotZ, xunit='ang', labPos=0.8*peak, lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-20)
  yr<-c(-peak,peak)
  hereFlux<-spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)]
  
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sky<-sky+abs(median(sky[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.45)
  
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sn<-sn+abs(median(sn[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.5)
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.4)
  tmp<-tmp+abs(median(tmp[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.3)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  
  range<-500
  wavePoint<-4150
  peak<-1.2*max(spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T)
  if (is.finite(peak)==F){peak=peak_F}
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=(1+plotZ)*wavePoint+c(-range,range), ylim=c(-peak,peak), main='Ca H&K, G-band', cex.axis=1.4)
  plotLines(z=plotZ, xunit='ang', labPos=0.8*peak, lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-20)
  yr<-c(-peak,peak)
  hereFlux<-spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)]
  
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sky<-sky+abs(median(sky[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.45)
  
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sn<-sn+abs(median(sn[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.5)
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.4)
  tmp<-tmp+abs(median(tmp[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.3)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  
  range<-250
  wavePoint<-5175
  peak<-1.2*max(spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T)
  if (is.finite(peak)==F){peak=peak_F}
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=(1+plotZ)*wavePoint+c(-range,range), ylim=c(-peak,peak), main='Mg', cex.axis=1.4)
  plotLines(z=plotZ, xunit='ang', labPos=0.8*peak, lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-20)
  yr<-c(-peak,peak)
  hereFlux<-spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)]
  
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sky<-sky+abs(median(sky[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.45)
  
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sn<-sn+abs(median(sn[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.5)
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.4)
  tmp<-tmp+abs(median(tmp[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.3)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  range<-250
  wavePoint<-5895
  peak<-1.2*max(spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T)
  if (is.finite(peak)==F){peak=peak_F}
  magplot(spec$wave, hanning.smooth(spec$flux, degree=degSmooth), xlab='Wavelength, Ang', ylab='Counts', grid=T, type='l', xlim=(1+plotZ)*wavePoint+c(-range,range), ylim=c(-peak,peak), main='Na', cex.axis=1.4)
  plotLines(z=plotZ, xunit='ang', labPos=0.8*peak, lty=2, cex=1.2, EmCol='blue', AbsCol='red', labOff=-20)
  yr<-c(-peak,peak)
  hereFlux<-spec$flux[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)]
  
  sky<-((spec$sky/(max(spec$sky,na.rm=T)-min(spec$sky,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sky<-sky+abs(median(sky[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.45)
  
  sn<-((spec$sn/(max(spec$sn,na.rm=T)-min(spec$sn,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.5)
  sn<-sn+abs(median(sn[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.5)
  
  tmp <- ((spec$tmpFlux/(max(spec$tmpFlux,na.rm=T)-min(spec$tmpFlux,na.rm=T)))*(max(hereFlux,na.rm=T)-min(hereFlux,na.rm=T))*0.4)
  tmp<-tmp+abs(median(tmp[which(spec$wave > ((1+plotZ)*wavePoint)-range & spec$wave < ((1+plotZ)*wavePoint)+range)],na.rm=T))-((yr[2]-yr[1])*0.3)
  
  lines(spec$tmpWave*(1+spec$z), tmp, col='indianred2')
  lines(spec$wave, sky, col='darkgreen')
  lines(spec$wave, sn, col='gold')
  
  par(mfrow = c(1, 1))
  par(mar=c(3.1,3.1,3.1,3.1))
  
  layout(matrix(c(1), 1, 1, byrow = TRUE))

options(warn=0)

}