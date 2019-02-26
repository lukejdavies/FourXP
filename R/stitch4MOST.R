#' Stitches together 4MOST arms 
#'
#' @description Function takes outputs formthe 4MOST ETC and filter response curves and stiches them togther. 
#' Can also scale to a desired magnitude in a given band. However, to do this you will need FourXP::getfilt installed
#' and have unpacked the calibrator files:
#' 
#' @description > LibPaths<-.libPaths()
#' @description > system(paste('tar -xvf ', LibPaths, '/FourXP/data/calibrators.tar --directory ', LibPaths, '/FourXP/data/',sep='')) 
#' 
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @param blue_file input FITS file for blue arm (output from 4MOST ETC)
#' @param green_file input FITS file for green arm (output from 4MOST ETC)
#' @param red_file input FITS file for red arm (output from 4MOST ETC)
#' @param blue_thru_file blue thougput file for ETC
#' @param green_thru_file green thougput file for ETC
#' @param red_thru_file red thougput file for ETC
#' @param filt_scale if you wish to scale to a magnitude in a particular band, name the band here as a string. Options can be found with ?getfilt
#' @param mag_scale AB magnitude in filt_scale to scale to
#' @param makePlot TRUE/FALSE make a PDF figure with the output spectrum and diagnositics
#' @param outName if makePlot==TRUE, string of output plot file names
#' @return A list containing wave, flux, error, scaled flux, scaled error, mag_scale, filt_scale for spliced spectrum
#' @examples 
#' blue_file<-"specout_template_TempSim_SF_Vhighz_LRS_blue.fits"
#' green_file<-"specout_template_TempSim_SF_Vhighz_LRS_green.fits"
#' red_file<-"specout_template_TempSim_SF_Vhighz_LRS_red.fits"
#' blue_thru_file<-"4FS-ETC_app/4FS_ETC_system_model_v0.2/LRS/lrs_blue_material_4fs_efficiency_total.fits"
#' green_thru_file<-"4FS-ETC_app/4FS_ETC_system_model_v0.2/LRS/lrs_green_material_4fs_efficiency_total.fits"
#' red_thru_file<-"4FS-ETC_app/4FS_ETC_system_model_v0.2/LRS/lrs_red_material_4fs_efficiency_total.fits"
#' filt_scale='VST_r'
#' mag_scale=19
#' makeplot=T
#' outName='stictch4MOSTout.pdf'
#' spec<-stitch4MOST(blue_file,green_file,red_file, blue_thru_file, green_thru_file, red_thru_file, filt_scale=filt_scale, mag_scale=mag_scale, makePlot=T, outName=outName)
#' @export
stitch4MOST<-function(blue_file,green_file,red_file, blue_thru_file, green_thru_file, red_thru_file, filt_scale='VST_r', mag_scale=NA, makePlot=T, outName='stictch4MOSTout.pdf'){
  
  
  # read throughput files
  blue_thru<-readFITS(blue_thru_file, hdu=1)
  green_thru<-readFITS(green_thru_file, hdu=1)
  red_thru<-readFITS(red_thru_file, hdu=1)
  
  # read data
  conditions<-readFITS(red_file, hdu=1)
  results_b<-readFITS(blue_file, hdu=2)
  results_g<-readFITS(green_file, hdu=2)
  results_r<-readFITS(red_file, hdu=2)
  
  # interpolate throuput to same wavelength scale as data
  red_thru2<-approx(red_thru$col[[1]]*10, red_thru$col[[12]], results_r$col[[1]])$y
  blue_thru2<-approx(blue_thru$col[[1]]*10, blue_thru$col[[12]], results_b$col[[1]])$y
  green_thru2<-approx(green_thru$col[[1]]*10, green_thru$col[[12]], results_g$col[[1]])$y
  
  # define band-specific fluxes
  red_flux<-(results_r$col[[3]])
  green_flux<-(results_g$col[[3]])
  blue_flux<-(results_b$col[[3]])
  
  # define band-specific wavelengths
  red_wave<-results_r$col[[1]]
  green_wave<-results_g$col[[1]]
  blue_wave<-results_b$col[[1]]
  
  # define band-specific nosie
  red_noise<-results_r$col[[4]]
  green_noise<-results_g$col[[4]]
  blue_noise<-results_b$col[[4]]
  
  # define band-specific S/N ratio
  red_SNR<-results_r$col[[5]]
  green_SNR<-results_g$col[[5]]
  blue_SNR<-results_b$col[[5]]
  
  # set up error vectors
  red_error<-red_flux
  green_error<-green_flux
  blue_error<-blue_flux
  
  
  ## add simulated error back in being randomly sampled from normal distribution with SD=Noise (suggested approch from Tom Dwelly)
  for (j in 1:length(blue_flux)){
    blue_error[j]<-rnorm(1,mean=0, sd=1.0*blue_noise[j])
    blue_flux[j]<-blue_flux[j]+blue_error[j]
  }
  for (j in 1:length(green_flux)){
    green_error[j]<-rnorm(1,mean=0, sd=1.0*green_noise[j])
    green_flux[j]<-green_flux[j]+green_error[j]
  }
  for (j in 1:length(red_flux)){
    red_error[j]<-rnorm(1,mean=0, sd=1.0*red_noise[j])
    red_flux[j]<-red_flux[j]+red_error[j]                                      
  }

  # divide through by throughput
  blue_fluxThru<- blue_flux/blue_thru2
  green_fluxThru<- green_flux/green_thru2
  red_fluxThru<-red_flux/red_thru2
  
  # set stupid values to zero
  blue_fluxThru[which(is.finite(blue_fluxThru)==F)]<-0
  green_fluxThru[which(is.finite(green_fluxThru)==F)]<-0
  red_fluxThru[which(is.finite(red_fluxThru)==F)]<-0
  

  #define common wavelength scale
  CommonWave<-seq(min(blue_wave,na.rm=T), max(red_wave,na.rm=T),0.3980254 )
  CommonFlux<-rep(NA,length(CommonWave))
  CommonError<-rep(NA,length(CommonWave))
  

  #rebin each arm to common wavelength scale
  blue_CommonFlux<-approx(blue_wave, blue_fluxThru, CommonWave)$y
  blue_CommonError<-approx(blue_wave, blue_error, CommonWave)$y
  green_CommonFlux<-approx(green_wave, green_fluxThru, CommonWave)$y
  green_CommonError<-approx(green_wave, green_error, CommonWave)$y
  red_CommonFlux<-approx(red_wave, red_fluxThru, CommonWave)$y
  red_CommonError<-approx(red_wave, red_error, CommonWave)$y
  
  

  # find overlapping regions in blue/green
  overLapBG<-CommonWave[which(CommonWave>min(green_wave,na.rm=T) &  CommonWave<=max(blue_wave,na.rm=T))]
  overLapBGPix<-which(CommonWave>min(green_wave,na.rm=T) &  CommonWave<=max(blue_wave,na.rm=T)) 
  scBG_B<-seq(1,0,length.out=floor(length(overLapBG)/2)) # define weighting in overlap region (1-0 for Blue)
  scBG_B<-c(rep(1,length(overLapBG)/4),scBG_B)
  tmp<-length(overLapBG)-length(scBG_B)
  scBG_B<-c(scBG_B, rep(0,tmp))
  scBG_G<-seq(0,1,length.out=floor(length(overLapBG)/2)) # define weighting in overlap region (0-1 for Green)
  scBG_G<-c(rep(0,length(overLapBG)/4),scBG_G)
  tmp<-length(overLapBG)-length(scBG_G)
  scBG_G<-c(scBG_G, rep(1,tmp))
  
  # find overlapping regions in green/red
  overLapGR<-CommonWave[which(CommonWave>min(red_wave,na.rm=T) &  CommonWave<=max(green_wave,na.rm=T))]
  overLapGRPix<-which(CommonWave>min(red_wave,na.rm=T) &  CommonWave<=max(green_wave,na.rm=T))
  scGR_G<-seq(1,0,length.out=floor(length(overLapGR)/2)) # define weighting in overlap region (1-0 for Green)
  scGR_G<-c(rep(1,length(overLapGR)/4),scGR_G)
  tmp<-length(overLapGR)-length(scGR_G)
  scGR_G<-c(scGR_G, rep(0,tmp))
  scGR_R<-seq(0,1,length.out=floor(length(overLapGR)/2)) # define weighting in overlap region (0-1 for Red)
  scGR_R<-c(rep(0,length(overLapGR)/4),scGR_R)
  tmp<-length(overLapGR)-length(scGR_R)
  scGR_R<-c(scGR_R, rep(1,tmp))
  
  # Scale both arms to meadin flux of green arem in overlap region
  scBlue<-median(blue_CommonFlux[overLapBGPix], na.rm=T)/median(green_CommonFlux[overLapBGPix], na.rm=T)
  scRed<-median(green_CommonFlux[overLapGRPix], na.rm=T)/median(red_CommonFlux[overLapGRPix], na.rm=T)
  blue_CommonFluxSc<-blue_CommonFlux/scBlue
  green_CommonFluxSc<-green_CommonFlux
  red_CommonFluxSc<-red_CommonFlux*scRed
  

  # Write relavent parts to CommonFlux vector for non-overalp regions.... 
  CommonFlux[which(CommonWave<min(green_wave,na.rm=T))]<-blue_CommonFluxSc[which(CommonWave<min(green_wave,na.rm=T))]
  CommonFlux[which(CommonWave>max(blue_wave,na.rm=T) & CommonWave<min(red_wave,na.rm=T))]<-green_CommonFluxSc[which(CommonWave>max(blue_wave,na.rm=T) & CommonWave<min(red_wave,na.rm=T))]
  CommonFlux[which(CommonWave>max(green_wave,na.rm=T))]<-red_CommonFluxSc[which(CommonWave>max(green_wave,na.rm=T))]
  # ...and then overlap region 
  CommonFlux[overLapBGPix]<-((blue_CommonFluxSc[overLapBGPix]*scBG_B)+(green_CommonFluxSc[overLapBGPix]*scBG_G))
  CommonFlux[overLapGRPix]<-((green_CommonFluxSc[overLapGRPix]*scGR_G)+(red_CommonFluxSc[overLapGRPix]*scGR_R))
  
  # Write relavent parts to CommonError vector... 
  CommonError[which(CommonWave<min(green_wave,na.rm=T))]<-blue_CommonError[which(CommonWave<min(green_wave,na.rm=T))]
  CommonError[which(CommonWave>max(blue_wave,na.rm=T) & CommonWave<min(red_wave,na.rm=T))]<-green_CommonError[which(CommonWave>max(blue_wave,na.rm=T) & CommonWave<min(red_wave,na.rm=T))]
  CommonError[which(CommonWave>max(green_wave,na.rm=T))]<-red_CommonError[which(CommonWave>max(green_wave,na.rm=T))]
  
  CommonError[overLapBGPix]<-((blue_CommonError[overLapBGPix]*scBG_B)+(green_CommonError[overLapBGPix]*scBG_G))
  CommonError[overLapGRPix]<-((green_CommonError[overLapGRPix]*scGR_G)+(red_CommonError[overLapGRPix]*scGR_R))
  
  # Write relavent parts to CommonFilter (weights) vector... 
  CommonFilterBlue<-rep(NA,length(CommonWave))
  CommonFilterBlue[which(CommonWave<min(green_wave,na.rm=T))]<-1
  CommonFilterBlue[overLapBGPix]<-scBG_B
  
  CommonFilterGreen<-rep(NA,length(CommonWave))
  CommonFilterGreen[which(CommonWave>max(blue_wave,na.rm=T) & CommonWave<min(red_wave,na.rm=T))]<-1
  CommonFilterGreen[overLapBGPix]<-scBG_G
  CommonFilterGreen[overLapGRPix]<-scGR_G

  CommonFilterRed<-rep(NA,length(CommonWave))
  CommonFilterRed[which(CommonWave>max(green_wave,na.rm=T))]<-1
  CommonFilterRed[overLapGRPix]<-scGR_R
  
  CommonFluxSc<-NA
  CommonErrorSc<-NA

  # This section scales the flux of the spectrum to the desired inpout magnitude and band
    if (is.na(mag_scale)==F){
      filter=getfilt(filt_scale)[,2:3]
      filt_wave <- filter[,1]
      filt_trans <- filter[,2]
      
      interp<-approx(filt_wave, filt_trans, CommonWave)
      filt_trans_interp <- interp$y
      filt_trans_interp[is.na(filt_trans_interp)] = 0
      
      wavefac=1e-10
      c<-299792458
      
      fluxnu=(wavefac*CommonFlux*CommonWave^2)/c
      
      flux_conv <- fluxnu*filt_trans_interp
      
      temp <- fluxnu*filt_trans_interp*CommonWave
      temp2 <- filt_trans_interp*CommonWave
      flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
      
      flux_need <- 10.^(-0.4*(mag_scale+48.6))
      
      flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
      CommonFluxSc <- (flux_sc_Hz*(2.998e18))/(CommonWave^2)
      CommonErrSc<-(CommonError/CommonFlux)*CommonFluxSc
      
    }
  
  
  # This section makes plots if requested
  if (makePlot==T){
    
    if (is.na(mag_scale)==T){
        
      pdf(outName, width=10,height=20)
      
      
      par(mfrow = c(6, 3))
      par(oma=c(0,0,0,0))
      par(mar=c(3.1,3.1,3.1,3.1))
      
      layout(matrix(c(1,2,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8), 6, 3, byrow = TRUE))
    }
    
    if (is.na(mag_scale)==F){
      
      pdf(outName, width=10,height=23)
      
      
      par(mfrow = c(7, 3))
      par(oma=c(0,0,0,0))
      par(mar=c(3.1,3.1,3.1,3.1))
      
      layout(matrix(c(1,2,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8 ,9,9,9), 7, 3, byrow = TRUE))
    }
    
    magplot(blue_wave, blue_flux, col='blue', ylab='Counts', xlab='Wavelength, Ang', type='l', main='Blue Arm')
    magplot(green_wave, green_flux, col='darkgreen', ylab='Counts', xlab='Wavelength, Ang', type='l', main='Green Arm')
    magplot(red_wave, red_flux, col='red', ylab='Counts', xlab='Wavelength, Ang', type='l', main='Red Arm')
    
    magplot(blue_wave, blue_fluxThru, type='l', col='blue', xlim=c(3800,9500), ylim=c(0,quantile(c(blue_fluxThru,green_fluxThru,red_fluxThru),probs=0.999,na.rm=T)), main='Arms X filter response', ylab='Counts', xlab='Wavelength, Ang')
    lines(green_wave, green_fluxThru, type='l', col='darkgreen')
    lines(red_wave, red_fluxThru, type='l', col='red')

    magplot(CommonWave, blue_CommonFluxSc, type='l', col='blue', xlim=c(3800,9500), ylim=c(0,quantile(c(blue_fluxThru,green_fluxThru,red_fluxThru),probs=0.999,na.rm=T)), main='Scaled to green in overlap', ylab='Counts', xlab='Wavelength, Ang')
    lines(CommonWave, green_CommonFluxSc, type='l', col='darkgreen')
    lines(CommonWave, red_CommonFluxSc, type='l', col='red')
    
    magplot(CommonWave, CommonFilterBlue, col='blue',type='l', ylab='Weight', xlab='Wavelength, Ang', main='Weighting Filter')
    lines(CommonWave, CommonFilterGreen, col='darkgreen')
    lines(CommonWave, CommonFilterRed, col='red')
  
    
    magplot(CommonWave, CommonFlux, type='l', col='black', xlim=c(3800,9500), ylim=c(0,quantile(c(blue_fluxThru,green_fluxThru,red_fluxThru),probs=0.999,na.rm=T)), main='Final Spliced', ylab='Counts', xlab='Wavelength, Ang')
    magplot(CommonWave, CommonError, type='l', main='Error Spectrum', ylab='Counts', xlab='Wavelength, Ang')
    
    
    if (is.na(mag_scale)==F){magplot(CommonWave, CommonFluxSc, type='l', col='black', xlim=c(3800,9500),  main=paste('Flux Scaled to ',mag_scale,'magAB in ',filt_scale, sep=''), ylab='Flux, ergs/sec/cm^2/Ang', xlab='Wavelength, Ang')}
    
    dev.off()
  
  }
  
  # returns list with relevent parts
  spec<-list(wave=CommonWave, flux=CommonFlux, error=CommonError, fluxSc=CommonFluxSc, errorSc=CommonErrorSc, mag_scale=mag_scale, filt_scale=filt_scale)
  return(spec)
  
  }