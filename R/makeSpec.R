#' Make a simulated spectrum using a set of SDSS+AGN templates at Air wavelengths
#'
#' @description Function takes inputs of redshift, magnitude, band for magnitude, g-i colour, stellar mass, 
#' sfr and AGN-type and returns a list with the paramters for a simulated spectrum. If AGN=='F': g-i colour, stellar mass 
#' and SFR are matched to a set of PPXF fits to SDSS spectra in the GAMA regions, the best-matched source spectrum is used.
#' This is then fit to the BC03 continuum models to extend in wavelength range. Note that the final flux sacling is not 
#' defined by col, mass and SFR, only the spectral type. If AGN='B'/'N', col, mass and SFR are ignored and either a 
#' broad-line ('B') or narrow ('N') line tempalte is used. This is still scaled to the correct magnitude and redshift. 
#' @param id id to assign to spectrum   
#' @param z redshfit for spectrum you wish to simulate
#' @param mag magnitude to scale final spectrum to
#' @param band photometric band in which scaling is to be done. See ?getfilt for options 
#' @param col g-i colour to match to template spectrum
#' @param mass Log10(stellar mass) to match to templale spectrum  
#' @param sfr SFR to match to templale spectrum 
#' @param agn Use AGN template? 'N'=no AGN template, 'B'= broad-line template, 'N'= narrow-line template
#' @param specDir location of model database - link to this will have been provided by running install4XP()  
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 

#' @export
makeSpec<-function(id='tmp', z=0.456, mag=19.8, band='VST_r', col=0.4, mass=10.2,sfr=10, agn='F', specDir='/Users/luke/work/IWG8/FourXPmodels/'){
  
  
  load(paste(specDir,'Spectra/modelSpecTab.Rdata', sep=''))
  modelTab<-data.frame(modelTab)
  tmp<-list.files(path=paste(specDir, 'Spectra/',sep=''), pattern='*')
  modelTab<-modelTab[which(modelTab$MODNAME %in% tmp==T),]
  logsSFR<-log10(modelTab$SFR/(10.^modelTab$LOGMSTAR))
  HA_EW<-modelTab$HA_EW
  
  
  spm<-read.table(paste(specDir,'BC03/csp_Zeq1p0_chab.sed',sep=''), header=F)
  age<-unique(spm[,1])
  spm2<-read.table(paste(specDir,'BC03/csp_Zeq0p2_chab.sed',sep=''), header=F)
  age2<-unique(spm2[,1])
  
  tmpArr<-array(NA,dim=c(6252,2,(length(age)+length(age2))))
  
  count<-1
  
  for (i in 1:ceiling((length(tmpArr[1,1,])/2.0))){
    
    spm_wave<-spm[which(spm[,1]==age),2]    
    spm_flux<-spm[which(spm[,1]==age),3]
    spm_flux<-spm_flux[which(spm_wave>800 & spm_wave<10100)]
    spm_wave<-spm_wave[which(spm_wave>800 & spm_wave<10100)]
    
    tmpArr[,1,count]<-spm_wave
    tmpArr[,2,count]<-spm_flux
    
    count<-count+1
  }
  
  for (j in 1:ceiling((length(tmpArr[1,1,])/2.0))){
    spm_wave<-spm2[which(spm2[,1]==age2[j]),2]    
    spm_flux<-spm2[which(spm2[,1]==age2[j]),3]
    spm_flux<-spm_flux[which(spm_wave>800 & spm_wave<10100)]
    spm_wave<-spm_wave[which(spm_wave>800 & spm_wave<10100)]
    
    tmpArr[,1,count]<-spm_wave
    tmpArr[,2,count]<-spm_flux
    
    count<-count+1
  }
  
  if (agn=='F'){
    
  
      
      #### Match Sim to Mod ####
      dist_col<-abs(col-modelTab$COL)/sd(modelTab$COL, na.rm=T)
      dist_mass<-abs(mass-modelTab$LOGMSTAR)/sd(modelTab$LOGMSTAR[which(is.finite(modelTab$LOGMSTAR)==T)], na.rm=T) 
      dist_sfr<-abs(log10(sfr)-log10(modelTab$SFR))/sd(log10(modelTab$SFR), na.rm=T)
      
      if (sfr>0) {
        dist_all<-sqrt(dist_col^2+dist_mass^2+dist_sfr^2)
        match<-which(dist_all==min(dist_all[which(logsSFR > -10.0 & HA_EW>=0.1)], na.rm=T))
      }
      if (sfr==0) {
        dist_all<-sqrt(dist_col^2+dist_mass^2)
        match<-which(dist_all==min(dist_all[which(logsSFR < -10.5 & HA_EW<0.1)]))
      }
      ##########################
      
      
      ### Get correct model ####       
      spec<-get.spec(paste(specDir,'Spectra/',modelTab$MODNAME[match],sep=''), xunit='ang', yunit='ang', z=0)
      
    }
    
    if (agn!='F'){
      agn_reds<-c(0.1,0.3,0.5,0.7,1.1,1.3,1.55,1.85,2.25,2.75,3.5,4.5,6.0)
      agn_redsN<-c('010','030','050','070','110','130','155','185','225','275','350','450','600')
      z_match<-agn_redsN[which(abs(z-agn_reds)==min(abs(z-agn_reds)))]
      z_match<-z_match[1]
      if (agn=='B'){name<-paste(specDir,'Spectra/AGN_v1.0_z',z_match,'_mr165_type1.fits',sep='')}
      if (agn=='N'){name<-paste(specDir,'Spectra/AGN_v1.0_z',z_match,'_mr165_type2.fits',sep='')}       
      spec<-get.spec(name, xunit='ang', yunit='ang', z=agn_reds[which(agn_redsN==z_match)])
      spec$wave<-AirFromVacuum(spec$wave)
      spec$wave<-spec$wave/(1+spec$z)
      spec$wave<-spec$wave-1.
      spec$z<-0
    }
    
    ### Scale spec in z ###
    specNew<-spec
    specNew$wave<-spec$wave*(1+z)
    specNew$z<-z
    
    ### Interpolate back to sensible wavelength range ###
    specNew$flux<-approx(specNew$wave,specNew$flux,seq(3000,10000,0.1))$y
    specNew$wave<-seq(3000,10000,0.1)
    
    specNew$flux[which(is.finite(specNew$flux)==F)]<-0.0
    
    
    ### Scale to desired r-mag ###
    wavefac=1e-10
    c<-299792458
    
    if (z<0.6) {filter=getfilt('VST_r')[,2:3]}
    if (z>=0.6) {filter=getfilt('VST_i')[,2:3]}
    filt_wave <- filter[,1]
    filt_trans <- filter[,2]
    
    interp<-approx(filt_wave, filt_trans, specNew$wave)
    filt_trans_interp <- interp$y
    filt_trans_interp[is.na(filt_trans_interp)] = 0
    
    fluxnu=(wavefac*specNew$flux*specNew$wave^2)/c
    
    flux_conv <- fluxnu*filt_trans_interp
    
    temp <- fluxnu*filt_trans_interp*specNew$wave
    temp2 <- filt_trans_interp*specNew$wave
    flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
    
    flux_need <- 10.^(-0.4*(mag+48.6))
    
    flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
    flux_sc <- (flux_sc_Hz*(2.998e18))/(specNew$wave^2)
    
    
    
    ##########################################
    #flux_sc[which(flux_sc==0 & wave_new<6000)]<-median(flux_sc[which(flux_sc>0)[1:50]],na.rm=T)
    #flux_sc[which(flux_sc==0 & wave_new>6000)]<-median(flux_sc[which(flux_sc>0)[(length(which(flux_sc>0))-50):length(which(flux_sc>0))]],na.rm=T)
    
    chi_min<-1e9
    
    for (k in 200:length(tmpArr[1,1,])){
      spm_wave<-tmpArr[,1,k]
      spm_flux<-tmpArr[,2,k]
      spm_wave<-spm_wave*(1+z)
      
      if (z<0.6) {filter=getfilt('VST_r')[,2:3]}
      if (z>=0.6) {filter=getfilt('VST_i')[,2:3]}
      
      
      filt_wave <- filter[,1]
      filt_trans <- filter[,2]
      
      interp<-approx(filt_wave, filt_trans, spm_wave)
      filt_trans_interp <- interp$y
      filt_trans_interp[is.na(filt_trans_interp)] = 0
      
      fluxnu=(wavefac*spm_flux*spm_wave^2)/c
      
      flux_conv <- fluxnu*filt_trans_interp
      
      temp <- fluxnu*filt_trans_interp*spm_wave
      temp2 <- filt_trans_interp*spm_wave
      flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
      flux_need <- 10.^(-0.4*(mag+48.6))
      
      flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
      spm_flux_sc <- (flux_sc_Hz*(2.998e18))/(spm_wave^2)
      
      
      sel<-which(flux_sc>0)
      flux_sc_interp<-approx(spm_wave,spm_flux_sc,specNew$wave[sel])$y
      chi<-sum((flux_sc_interp-flux_sc[sel])^2/flux_sc[sel])
      
      if (chi<chi_min){
        sel2<-k
        chi_min<-chi
      }
      
      
    }
    
    specComp<-list(wave=specNew$wave,flux=flux_sc, xunit='ang', yunit='ang')
    i_mag<-magABspec(specComp,filter='VST_i')
    
    spm_wave<-tmpArr[,1,sel2]
    spm_flux<-tmpArr[,2,sel2]
    spm_wave<-spm_wave*(1+z)
    
    filter=getfilt('VST_i')[,2:3]
    filt_wave <- filter[,1]
    filt_trans <- filter[,2]
    
    interp<-approx(filt_wave, filt_trans, spm_wave)
    filt_trans_interp <- interp$y
    filt_trans_interp[is.na(filt_trans_interp)] = 0
    
    fluxnu=(wavefac*spm_flux*spm_wave^2)/c
    
    flux_conv <- fluxnu*filt_trans_interp
    
    temp <- fluxnu*filt_trans_interp*spm_wave
    temp2 <- filt_trans_interp*spm_wave
    flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
    flux_need <- 10.^(-0.4*(i_mag+48.6))
    
    flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
    spm_flux_sc <- (flux_sc_Hz*(2.998e18))/(spm_wave^2)
    
    
    if (length(which(flux_sc==0 & specNew$wave<specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])>0)) {
      lowCut<-max(specNew$wave[which(flux_sc==0 & specNew$wave<specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])])
      low_sel<-which(specNew$wave<=(lowCut+500))
      spm_low_sel<-which(spm_wave<=(lowCut+500))
      spm_low_interp<-approx(spm_wave,spm_flux_sc,specNew$wave[low_sel])$y
      off_low<-median(spm_low_interp[which(flux_sc[low_sel]>0)]-flux_sc[which(flux_sc[low_sel]>0)],na.rm=T)
      spm_low_interp<-spm_low_interp-off_low
      low_sel2<-which(specNew$wave<=(lowCut))
      flux_sc[low_sel2]<-spm_low_interp[low_sel2]
    }
    
    if (length(which(flux_sc==0 & specNew$wave>specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])>0)) {
      highCut<-min(specNew$wave[which(flux_sc==0 & specNew$wave>specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])])
      high_sel<-which(specNew$wave>=(highCut-500))
      spm_high_sel<-which(spm_wave>=(highCut-500))
      spm_high_interp<-approx(spm_wave,spm_flux_sc,specNew$wave[high_sel])$y
      off_high<-median(spm_high_interp[which(flux_sc[high_sel]>0)]-flux_sc[high_sel[which(flux_sc[high_sel]>0)]],na.rm=T)
      spm_high_interp<-spm_high_interp-off_high
      high_sel2<-which(specNew$wave>=(highCut))
      flux_sc[high_sel2]<-approx(specNew$wave[high_sel],spm_high_interp,specNew$wave[high_sel2])$y
    }
    
    
    
    ##########################################

    flux_sc<-flux_sc+abs(min(flux_sc, na.rm=T))

    
    wavefac=1e-10
    c<-299792458
    
    filter=getfilt(band)[,2:3]      
    filt_wave <- filter[,1]
    filt_trans <- filter[,2]
    
    interp<-approx(filt_wave, filt_trans, specNew$wave)
    filt_trans_interp <- interp$y
    filt_trans_interp[is.na(filt_trans_interp)] = 0
    
    fluxnu=(wavefac*flux_sc*specNew$wave^2)/c
    
    flux_conv <- fluxnu*filt_trans_interp
    
    temp <- fluxnu*filt_trans_interp*specNew$wave
    temp2 <- filt_trans_interp*specNew$wave
    flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
    flux_need <- 10.^(-0.4*(mag+48.6))
    
    flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
    flux_sc <- (flux_sc_Hz*(2.998e18))/(specNew$wave^2)
    
  spec<-list(id=id, wave=specNew$wave, flux=flux_sc, z=z, mag=mag, band=band, col=col, mass=mass,sfr=sfr,agn=agn, expMin=NA, prob=NA)
 return(spec)
  
  
}
