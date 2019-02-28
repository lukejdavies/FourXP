#' 4MOST redshifting tool for 1D spectra
#'
#' @description Based on redshfiting tool originaly developed by Ivan Baldry in IDL (Baldry et al., 2014 MNRAS). 
#' Repurposed & expanded in R by Luke Davies and Leon Drygala. Takes a 1D spectrum and perfoms a 
#' Fourier cross corellation with a template spectra set. Returns the best fit redshifts
#' and probabilities.
#' 
#' @param specRaw R struture contianing information required information.
#' Must have the following componenents: specRaw$wave=vector of spectrum wavelengths, 
#' specRaw$flux=vector spectrum fluxes the same size as specRaw$wave, specRaw$error=vector 
#' spectrum errors the same size as specRaw$wave, and if doHelio=T, specRaw$RA=spectrum RA,
#' specRaw$DEC=spectrum DEC, specRaw$UTMJD = observation Jullian date, specRaw$longitude = 
#' observatory longitude, specRaw$latitude = observatory latitude, specRaw$altitude = 
#' observatory altitude.    
#' @param tempFile Path to file containing spectral template data 
#' @param oversample wavelength oversampling rate
#' @param num number of crosscorrelation peaks to identify
#' @param templateNumbers template numbers to use in fitting
#' @param stLambda lower bound of the wavelength range to fit over 
#' @param endLambda = upper bound of the wavelength range to fit over
#' @param minval minmum value to reject croos correlations
#' @param maxval maximum value to reject croos correlations
#' @param z_prior redshift prior, two element vector with c(lo, hi)
#' @param doHelio TRUE/FALSE perform helocentric correction. If TRUE you must 
#' provide RA,DEC,UTMJD, longitude, latitude and altitue in the specRaw structure. 
#' @param verbose  TRUE/FLASE - let me know what's going on.
#' @return a list containing various outputs from the redshifting code. key parameters are: 
#' results$Z=best-fit redshift, results$Z1_PROB=probability of best-fit redshift and results$TEMPLATE=best fit template number. 
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plot(spec$wave, hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts')
#' tempFile=paste(.libPaths(),'/FourXP/data/calibrators/AutoZTemp/filtered-templates.fits',sep='')
#' FourXP_Z_out<-FourXP_Z(spec, tempFile=tempFile, doHelio=F)
#' plotLines(z=FourXP_Z_out$z)
#' cat('Probability of correct redshift is: ', FourXP_Z_out$prob, '\n')
#' @export
FourXP_Z= function(specRaw, tempFile = NA,oversample = 5, num = 5, templateNumbers = c(2:14,16:22,40:47), stLambda = 3726, endLambda = 8850, minval = -1.0e4, maxval = 1.0e6, z_prior=c(-1,1000), doHelio=T,highZ=T, verbose = TRUE){
  
  if (is.null(specRaw$error)){
    
    cat('**WARNING** No error supplied! Using dummy error....','\n')
    specRaw$error<-rep(median(specRaw$flux, na.rm=T)/10000, length(specRaw$flux))
    
  }
  
  if (is.na(tempFile)){tempFile = paste(.libPaths(),'/FourXP/data/calibrators/AutoZTemp/filtered-templates.fits',sep='')}
  
  # TODO remove timing
  TOTALTIME <- proc.time()
  lambda<-specRaw$wave
  useInvCorrection = TRUE
  specNum<-1
  
  #set up new lambda scale to rebin spectrum and templates to
  logLambdaData <- SetUpLogLambda(verbose = verbose, oversample = oversample, highZ=highZ)
  newLogLambda <- logLambdaData$logLambda
  
  
    specRaw$lambda <- lambda
PROCESSTIME <- proc.time()
  spec <- ProcessSpectrum(specRaw, stLambda = stLambda, endLambda = endLambda, minval = minval, maxval = maxval, 
                          useInvCorrection = useInvCorrection,  verbose = verbose)
  PROCESSTIME <- proc.time()[3] - PROCESSTIME[3]
  
  spec$countHighval  <- length(which(spec$flux > 20))
  spec$countLowval   <- length(which(spec$flux < -20))
  spec$meanadNorm    <- mean(abs(spec$flux))
  spec$rmsNorm       <- sqrt( mean(spec$flux^2) )
  
  
  # Convert spectral lambda to vacuum wavelength from air wavelength input.
  logVlambda <- log(VacuumFromAir(spec$lambda),10)  
  length <- length(logVlambda)
  # Rebin filtered spectrum ensuring zero outside range.
  if (highZ==F) {specRebin <- approx(x = c(3.3, logVlambda[1] - 0.001, logVlambda, logVlambda[length]+0.001, 4.0), y = c(0, 0, spec$flux, 0, 0), xout = newLogLambda, method = "linear", yleft=0.0, yright=0.0)}#TODO check interpol is correct
  if (highZ==T) {specRebin <- approx(x = c(3.0, logVlambda[1] - 0.001, logVlambda, logVlambda[length]+0.001, 4.0), y = c(0, 0, spec$flux, 0, 0), xout = newLogLambda, method = "linear", yleft=0.0, yright=0.0)} #TODO check interpol is correct
  spec$lambda <- specRebin$x
  spec$flux <- specRebin$y
  
  # load rebinned template data 
  tempData <- RebinTempData(newLogLambda, templateNumbers=templateNumbers, 
                                                   file = tempFile, verbose = verbose)
  
  helioVel<-0
  
  if (doHelio==T) {helioVel <- Heliocentric(spec$RA*180/pi, spec$DEC*180/pi, 2000, jd = spec$UTMJD, longitude = spec$longitude, 
                          latitude = spec$latitude, altitude = spec$altitude)}else{helioVel<-0}

 plan = 0
CROSSTIME <- proc.time()
  #get cross correlation info and find highest peaks in data
  ccinfo <- DoCrossCorr(spec = spec, gap = logLambdaData$gap, tempData = tempData, helioVel = helioVel, plan = plan, z_prior=z_prior,highZ=highZ)
CROSSTIME <- CROSSTIME[3] - proc.time()[3]
  peaks <- FindHighestPeaks(ccinfo, num=num)
  #fit quadratic and adjust redshift slightly. Also save ccSigma data
  spec$ccSigma <- rep(-1,length(peaks))
  i <- 0
  for( peak in peaks ){
      
    peak$redshift   <- FindMax(xArray = ccinfo[[peak$templateID]]$shifts[(peak$shiftIndex-3):(peak$shiftIndex+3)],
                               yArray = ccinfo[[peak$templateID]]$crossCorrRaw[(peak$shiftIndex-3):(peak$shiftIndex+3)],
                               n = 2)
    spec$ccSigma[i <- i+1] <- peak$crossCorr
  }
  spec$z <- peaks[[1]]$redshift
  #calculate FOM and Probability
  spec <- CalculateCertainty(spec)
  
  #print results
  if(verbose){
    cat('\nRedshifts for spectrum number', specNum,'are as follows\n')
    cat("Template\tRedshift\t\tCrossCorr\t\tShiftIndex\n")
    for( peak in peaks ){
      cat(peak$template, peak$redshift, peak$crossCorr, peak$shiftIndex,"\n", sep="\t\t")
    }
    cat("Probabilty on best match is ",spec$prob,"\n")
  }
  
  spec$results        <- c(peaks[[1]]$redshift, spec$prob, peaks[[1]]$crossCorr, peaks[[1]]$template, peaks[[2]]$redshift, 
                           peaks[[2]]$crossCorr, peaks[[2]]$template, peaks[[3]]$crossCorr, peaks[[4]]$crossCorr)
  names(spec$results) <- c('Z', 'Z1_PROB', 'CC_SIGMA','TEMPLATE','Z2','CC_SIGMA2','TEMPLATE2','CC_SIGMA3','CC_SIGMA4')
  TOTALTIME <- proc.time()[3] - TOTALTIME[3]
  spec$timings <- c(TOTALTIME, CROSSTIME, PROCESSTIME)
  names(spec$timings) <- c('totalTime', 'crossTime', 'proccessTime')
  return = spec
}






