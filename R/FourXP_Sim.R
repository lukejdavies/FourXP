#' Simulated and mock observe spectra with 4MOST
#'
#' @description High level function that combines most of the functionality of FourXP. 
#' Allows the user to generate a simulated spectrum and 'observe' with the 4MOST ETC. 
#' Various statisitics about the observation are returned. 
#' 
#' @param id ID to assing to this simulation
#' @param zIn input redshift to simulate 
#' @param mag magnitude to scale final spectrum to
#' @param band photometric band in which scaling is to be done. See ?getfilt for options 
#' @param col g-i colour to match to template spectrum
#' @param mass Log10(stellar mass) to match to templale spectrum  
#' @param sfr SFR to match to templale spectrum 
#' @param agn Use AGN template? 'F'=no AGN template, 'B'= broad-line template, 'N'= narrow-line template
#' @param specDir location of model database, this is required but is somewhat large (~2.5Gb). As such,it can be found here: https://www.dropbox.com/s/hwnbabirn59796x/FourXPmodels.zip?dl=0    
#' @param expMin Total exposure time to simulate observation  
#' @param nSub number of sub-exposure in expMin
#' @param SKYBRIGHT_TYPE where to measuresky brightness ZENITH or LOCAL?
#' @param AIRMASS simulated airmass of observation
#' @param IQ image quality value to simulate (V-band,FWHM,arcsec) 
#' @param SKYBRIGHT sky brightness in mag/arcsec^2
#' @param TILT fibre tilt in mm
#' @param MISALIGNMENT fibre misallignment in arcsec
#' @param systemModDir IMPORTANT: directory location of 4FS system files. Will be something like '4FS-ETC_app/4FS_ETC_system_model_v0.2/' 
#' @param plot TRUE/FALSE plot the output spectrum with plotSpec()
#' @param verbose TRUE/FALSE let me know what's going on 
#' @return a list contining lots and lots of things from the simulated spectrum 
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' out<-FourXP_Sim(id='Test', zIn=0.234, mag=19.8, band='VST_r', col=0.55, mass=10.2,sfr=7, agn='F', expMin=60, nSub=3, verbose=TRUE)
#' @export
FourXP_Sim<-function(id='Test', zIn=0.234, mag=19.8, band='VST_r', col=0.55, mass=10.2,sfr=7,agn='F', specDir='/Users/luke/work/IWG8/FourXPmodels/', expMin=60, nSub=3, SKYBRIGHT_TYPE='ZENITH', AIRMASS=1.4, IQ=1.1, SKYBRIGHT=21.77, TILT=6.0, MISALIGNMENT=0.1, systemModDir='/Applications/4FS-ETC_app/4FS_ETC_system_model_v0.2/', plot=T, verbose=TRUE){
  
  cat('Running FourXP_Sim, please wait.....', '\n')
  if (verbose==T){cat('     - Making simulated spectrum...', '\n')}
  spec<-makeSpec(id=id, z=zIn, mag=mag, band=band, col=col, mass=mass,sfr=sfr,agn=agn, specDir=specDir)
  if (verbose==T){cat('     - Observing spectrum using 4FS ETC...', '\n')}
  observeSpec4FSOut<-observeSpec4FS(spec, expMin=expMin, nSub=nSub, keepFITS=F, SKYBRIGHT_TYPE=SKYBRIGHT_TYPE, AIRMASS=AIRMASS, IQ=IQ, SKYBRIGHT=SKYBRIGHT, TILT=TILT, MISALIGNMENT=MISALIGNMENT, systemModDir=systemModDir)
  if (verbose==T){cat('     - Stictching Spectral Arms...', '\n')}
  specObs<-stitch4MOST(observeSpec4FSOut=observeSpec4FSOut)
  if (verbose==T){cat('     - Running 4XP_Z...', '\n')}
  FourXP_ZOut<-FourXP_Z(specObs, verbose=F, doHelio=F)
  
  
  
  specObs$id<-id
  specObs$z<-FourXP_ZOut$z
  res<-(2.998e5)*abs(specObs$z-zIn)/(1+zIn)
  specObs$zIn<-zIn
  specObs$prob<-FourXP_ZOut$prob
  specObs$res<-res
  specObs$col<-col
  specObs$mass=mass
  specObs$sfr=sfr
  specObs$agn=agn
  
  if (verbose==T){
    
    cat('\n', 'KEY INPUT PROPERTIES:','\n \n')
    cat('   ID = ',id,'\n')
    cat('   Input Redshift = ',zIn,'\n')
    cat('   ABMag = ',mag,'\n')
    cat('   MagBand = ',band,'\n')
    if (agn=='F') {cat('   Template g-i col = ',col,'\n')}
    if (agn=='F') {cat('   Template Log[M*] = ',mass,'\n')}
    if (agn=='F') {cat('   Template SFR = ',sfr,'\n')}
    if (agn=='B') {cat('   AGN Type = Broad-line','\n')}
    if (agn=='B') {cat('   AGN Type = Narrow-line','\n')}
    cat('   Explorue Time (min) = ',expMin,'\n')
    cat('   Number of sub-exposures = ',nSub,'\n')
    cat('   Airmass = ',AIRMASS,'\n')
    cat('   Sky Brightness = ',SKYBRIGHT,'\n')
    
    cat('\n', 'KEY OUTPUT PROPERTIES:','\n \n')
    
    cat('Measured Redshift = ', specObs$z, '\n')
    cat('Measured Probability = ', specObs$prob, '\n')
    cat('Redshift Precision (in vs measured) = ', res, 'km/s \n')
    cat('Signal to Noise Blue (median 4200-5000) = ', median(specObs$blueRawSNR[which(specObs$blueRawWave>4200 & specObs$blueRawWave<5000)], na.rm=T), '\n')
    cat('Signal to Noise Green (median 5800-6600) = ', median(specObs$greenRawSNR[which(specObs$greenRawWave>5800 & specObs$greenRawWave<6600)], na.rm=T), '\n')
    cat('Signal to Noise GRed (median 7800-8600) = ', median(specObs$redRawSNR[which(specObs$redRawWave>7800 & specObs$redRawWave<8600)], na.rm=T), '\n')
    
    cat('\n \n')
  }
    
  cat('FourXP_Sim finished!', '\n')
  
  if (plot==T){
    specObsTmp<-specObs
    specObsTmp$flux<-specObs$fluxSc
    plotSpec(specObs)
  }
  
  return(specObs)
  
  
  
  
}