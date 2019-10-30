scFlux<-function(spec, mag, band){

wavefac=1e-10
    c<-299792458
    
    filter=getfilt(band)[,2:3]      
    filt_wave <- filter[,1]
    filt_trans <- filter[,2]
    
    interp<-approx(filt_wave, filt_trans, spec$wave)
    filt_trans_interp <- interp$y
    filt_trans_interp[is.na(filt_trans_interp)] = 0
    
    fluxnu=(wavefac*spec$flux*spec$wave^2)/c
    
    flux_conv <- fluxnu*filt_trans_interp
    
    temp <- fluxnu*filt_trans_interp*spec$wave
    temp2 <- filt_trans_interp*spec$wave
    flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
    flux_need <- 10.^(-0.4*(mag+48.6))
    
    flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
    flux_sc <- (flux_sc_Hz*(2.998e18))/(spec$wave^2)

    spec$flux<-flux_sc
    spec$mag<-mag
    
    return(spec)
    
}
