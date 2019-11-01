#' Processes spectrum to correct FourXP format
#'
#' @description Process spectrum to remove bad data, etc
#' 
#' @param specIn input spectrum
#' @param stLambda min wavelength range to fit over
#' @param endLambda max wavelength range to fit over
#' @param os1 length to extend below stLambda
#' @param os2 length to extend above endLambda
#' @param minval min value to reject y values
#' @param maxval max value to reject y values
#' @param clipvalue value to itterative clip peaks
#' @param baddataErrorValue value to replace bad data with
#' @param useInvCorrection set relative flux calibration for GAMA 
#' @param medWidth pixel  width of filtering kernal for FitandFilter -> MedSmooth 
#' @param smoothWidth pixel width of smoothing kernal for FitandFilter -> MedSmooth
#' @param runningWidth pixel width of smoothing kernal for FitandFilter ->RunningMeanSmooth
#' @param verbose TRUE/FALSE,  tell me whats going on
#' @author I. Baldry, L. Davies
#' @examples 
#' None applicable
#' @export
ProcessSpectrum = function(specIn, stLambda = 3726, endLambda = 8850, os1 = 10, os2 = 60, 
                           minval = -1.0e4, maxval = 1.0e6, clipvalue = 25, baddataErrorValue = 1e10, 
                           useInvCorrection = TRUE, medWidth=51, smoothWidth= 121, runningWidth=21, verbose = TRUE){
  
  
  spec <- specIn
  # Criteria for bad spectral points: non-finite - flagged as bad - 
  #   or large deviations of sky-subtracted data - assumed bad pixels.
  criteriaBadspec <- !is.finite(spec$flux) | (spec$flux < minval) | (spec$flux > maxval)
  # Criteria for bad errors: non-finite or <= zero
  criteriaBaderror = !is.finite(spec$error) | (spec$error <= 0)
  
  # Set approximate relative flux 'calibration' as a function of microns.
  # Was obtained usng GAMA standard stars. 
  if (useInvCorrection) {
    res <- c(  -3.406,  22.175, -49.208,  55.037, -23.441)
    invCorrection <- PolyCalc(res, spec$lambda/10000)
    # Apply relative flux calibration.
    spec$flux <- spec$flux / invCorrection
    spec$error <- spec$error / invCorrection
  }
  
  # Robustly high-pass filter the spectrum.
  r <- FitandFilter(spec, stLambda, endLambda, os1, os2, verbose = verbose)
  spec <- r[[1]]
  useFit <- r[[2]]
  
  
  useOkError <- which(!criteriaBaderror)
  # Broaden error over all lines
  spec$error[useOkError] <- AdjustBroaden(spec$error[useOkError])
  # Set minimum error to avoid anomalously low errors, i.e. take the larger of adjustedError and spec$error
  adjustedError <- 0.7*MedianAdjust(spec$error[useOkError], 13)
  booArray <- adjustedError > spec$error[useOkError]
  spec$error[useOkError][booArray] <- adjustedError[booArray]
  
  
  
  # Set bad-data error value where there exists any bad data .
  useBaddata <- which(criteriaBadspec | criteriaBaderror)
  spec$countBadPix <- length(useBaddata)
  spec$error[useBaddata] <- baddataErrorValue
  spec$flux[useBaddata] <- 0.0
  
  
  
  # divide by sigma^2 as suggested by Saunders, Cannon and Sutherland
  spec$flux = spec$flux / spec$error^2
  # clip value, was 30. 
  spec$flux = NormaliseMeandev(spec$flux, clipvalue=clipvalue, use=useFit)
  
  return = spec
}
