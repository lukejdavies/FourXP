#' Internal AutoZ function for setting up log wavelength scale
#'
#' @description Set up log wavelength scale
#' 
#' @param oversample how much to oversample the spectrum given the wavelength gap 'sdssGap'
#' @param lowValue min values of scale
#' @param highValue max value of scale
#' @param verbose TRUE/FALSE, tell me what's going on
#' @return logLambda wavelength scale
#' @author I. Baldry
#' @examples 
#' LamScale<-SetUpLogLambda(oversample=5, lowValue=3.3, highValue=4.0, verbose = TRUE)
#' @export
SetUpLogLambda = function(oversample=5, lowValue=3.3, highValue=4.0, verbose = TRUE,highZ=T){
  
  if (highZ==T){
    lowValue<-3.0
    highValue<-4.0
  }
  
  sdssGap <- 0.0001
  gap <- NULL
  # oversample the SDSS log wavelength gap, default is a factor of 5
  if (oversample != 0){
    gap <- sdssGap / oversample
  } #TODO throw div 0 exception
  
  velocityGap <- gap*log(10)*2.998E5
  if(verbose) cat('\nSET_UP_logLambda: Velocity rebinned gap', velocityGap)
  num <- round( (highValue - lowValue)/gap )

  logLambda <- lowValue + gap*(0:(num-1))
  
  if(verbose){ 
    cat('\nSET_UP_logLambda: Num points and range', num, 
      10.^logLambda[0], 10.^logLambda[num-1])
  }
  
  return = list("logLambda" = logLambda, "gap" = gap, "velocity_gap" = velocityGap)
}
