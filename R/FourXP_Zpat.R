#' Test version of redshfiting code for emission line pattern matching 
#'
#' @description Identifies emission lines in target spectrum and 
#' pattern matches to a catlogue of known emission line. 
#' 
#' @param specRaw Input spectrum, list in format of all of the things in FourXP 
#' @param z_prior redshift prior
#' @param plotCorr TRUE/FALSE plot the output of the match
#' @param filter TRUE/FALSE filter the spectrum for noise line peaks
#' @param verbose TRUE/FALSE let me know what's going on 
#' @return a list contining lots and lots of things from the simulated spectrum 
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' spec<-makeSpec(id='D134', z=zIn, mag=21, band='VST_r', col=0.6, mass=11,sfr=10, agn='F')
#' zFit<-FourXP_Zpat(spec, z_prior=c(zIn-0.1,zIn+0.1), plotCorr=T)
#' @export
FourXP_Zpat<-function(specRaw, z_prior=c(-1,5), plotCorr=T, verbose = TRUE, filter=TRUE){
  

  splineFit<-smooth.spline(specRaw$wave,specRaw$flux, df=7)
  spec<-specRaw
  spec$flux<-spec$flux-approx(splineFit$x, splineFit$y, spec$wave)$y
  spec2<-spec
  if (filter==T){
    for (j in 1:2){
      for (i in 3:length(spec$flux)-2){
       if (spec$flux[i]>1.5*sd(spec$flux, na.rm=T) & min(c(spec$flux[(i-2):(i-1)],spec$flux[(i+1):(i+2)]),na.rm=T)<0.5*sd(spec$flux, na.rm=T)){spec2$flux[i]=0}
       if ((-1*spec$flux[i])>1.5*sd(spec$flux, na.rm=T) & min(-1*c(spec$flux[(i-2):(i-1)],spec$flux[(i+1):(i+2)]),na.rm=T)<0.5*sd(spec$flux, na.rm=T)){spec2$flux[i]=0}
      }
    spec<-spec2
    }
  }

  pixArr<-seq(1,length(spec$flux),1)
  tmp<-spec$flux[which(spec$flux<2.0*sd(spec$flux,na.rm=T))]
  
  peaks<-which(spec$flux>2.5*sd(tmp,na.rm=T))
  

  peaks<-peaks[which(spec$wave[peaks]<7550 | spec$wave[peaks]>7750)]
  peaks<-peaks[which(peaks>100 & peaks<(length(spec$wave)-100))]
  
  tmpArr<-rep(0,length(spec$flux))
  tmpArr[peaks]<-1
  
  on<-0
  peakNew<-c()
  for (i in 3:(length(tmpArr)-2)){
    if (on==0 & sum(tmpArr[(i-3):(i-1)],na.rm=T)==0 & sum(tmpArr[(i):(i+2)],na.rm=T)>0){
      on=1
      onpos=i
      }
    if (on==1 & sum(tmpArr[(i-3):(i-1)],na.rm=T)>0 & sum(tmpArr[(i):(i+2)],na.rm=T)==0){
      on=0
      peakNew<-c(peakNew, ((i-1)+onpos)/2)
      onpos<-NA
      }
  }
 
  peakNew2<-peakNew
  for (i in 1:length(peakNew)){
    tmp<-abs(spec$flux[(peakNew[i]-10):(peakNew[i]+10)])
    tmp2<-pixArr[(peakNew[i]-10):(peakNew[i]+10)]
    sel<-which(tmp==max(tmp,na.rm=T))
    peakNew2[i]<-pixArr[tmp2[sel]]
  }
  
 
  peaks<-spec$wave[peakNew2]

  
  
  
  peakStrength<-abs(spec$flux[peakNew])
  peakStrength<-peakStrength/max(peakStrength,na.rm=T)
  
  lines<-load.lines()
  lines$wave_ang<-AirFromVacuum(lines$wave_ang)
  linesAng<-lines$wave_ang
  zTry<-seq(z_prior[1],z_prior[2],0.0001)
  offAll<-rep(NA, length(zTry))
  for (i in 1:length(zTry)){
    off<-rep(NA, length(linesAng))
    for (j in 1:length(linesAng)){
      linesZtry<-linesAng[j]*(1+zTry[i])
      if (min(abs(linesZtry-peaks),na.rm=T)<20){
        off[j]<-min(abs(linesZtry-peaks),na.rm=T)*(1.1-peakStrength[which(abs(linesZtry-peaks)==min(abs(linesZtry-peaks),na.rm=T))])*20*sqrt((1.1-(lines$strength[j]/max(lines$strength))))
      } else {
          off[j]<-20
          }
    }
    
    
    offAll[i]<-sum(off, na.rm=T)/length(which(off!=20))
    
  }
  
  zMatch<-zTry[which(offAll==min(offAll, na.rm=T))]
  peakBest <- -1e6
  for (i in 1:length(peaks)){
    if (min(abs((linesAng*(1+zMatch))-peaks[i]), na.rm=T)<20){
      if (spec$flux[peaks[i]]>peakBest){
        lineUse<-which(abs((linesAng*(1+zMatch))-peaks[i])==min(abs((linesAng*(1+zMatch))-peaks[i]), na.rm=T))
        peakUse<-peaks[i]
        peakBest<-spec$flux[peaks[i]]
        }
    }
  }
  
  zMatch<-(peakUse/linesAng[lineUse])-1
  

  
  if (plotCorr==T){
    
    par(mfrow = c(2, 1))
    par(mar=c(3.1,3.1,1.1,1.1))
    
    layout(matrix(c(1,2), 2, 1, byrow = TRUE))
    
    plotSpec.basic(spec, degSmooth=7)
    plotLines(z=specRaw$z, Air=T)
    points(peaks, rep(max(spec$flux,na.rm=T)*0.9,length(peaks)), col='gold', pch=16)
    points(linesAng*(1+zMatch), rep(max(spec$flux,na.rm=T)*0.85,length(linesAng)), col='darkgreen', pch=16)
    legend('bottomright', legend=c('Peaks Used', 'Line-list at fitted Z'), text.col=c('gold', 'darkgreen'))

    magplot(zTry, offAll, xlab='Redshfit',ylab='Correlation', type='l', xlim=c(zMatch-0.1, zMatch+0.1))
    lines(c(zMatch,zMatch), c(0, 1e10), col='red', lwd=2, lty=1)
    lines(c(specRaw$z,specRaw$z), c(0, 1e10), col='blue', lwd=2, lty=2)
    legend('bottomright', legend=c('zIn', 'zFit'), text.col=c('blue', 'red'))
  }
  
  if (verbose==T){
    res<-(2.998e5)*abs(zMatch-specRaw$z)/(1+specRaw$z)
    cat('Measured Redshift = ', zMatch, '\n')
    cat('Redshift Precision (in vs measured) = ', res, 'km/s \n')
  }
  
    
  return(zMatch)
  
}
