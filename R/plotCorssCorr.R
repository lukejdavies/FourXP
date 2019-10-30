#' Plot the cross-correlation output from FourXP_Z()
#'
#' @description Takes the output from FourXP_Z() and allows the user to plot the cross-correlation for a given template 
#' 
#' @param spec R struture output from FourXP_Z() 
#' @param templateNum Template number to plot correlation for 
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plot(spec$wave, hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts')
#' tempFile=paste(.libPaths(),'/FourXP/data/calibrators/AutoZTemp/filtered-templates.fits',sep='')
#' FourXP_Z_out<-FourXP_Z(spec, tempFile=tempFile, doHelio=F)
#' plotCrossCorr(FourXP_Z_out, templateNum=FourXP_Z_out$results[4])
#' @export 
plotCrossCorr<-function(spec, templateNum=13){
  
  ccinfo<-spec$ccinfo
  lengthTemp<-dim(summary(ccinfo))[1]
  numTemp<-rep(NA,lengthTemp)
  for (i in 1:lengthTemp){
    numTemp[i]<-as.numeric(ccinfo[[i]][1])
  }
  sel<-which(numTemp==templateNum)
  if (length(sel)==0){
    cat('**** WARNING: template number is not valid, please check you ran FourXP_Z() with this template ****', '\n')
    return(NULL)
    
  }
  magplot(ccinfo[[sel]]$shifts, ccinfo[[sel]]$crossCorr, type='l', xlab='Redshift', ylab='Cross Correlation')
  lines(c(spec$results[1], spec$results[1]), c(-500,500), col='red', lty=2)
  legend('topright', legend='Best Redshift', lty=2, col='red')
  
}
