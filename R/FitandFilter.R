#' FitandFilter
#' @description Set up region to use for fitting between start Lambda and end Lambda.
#' @param spec list with elements spec$flux and spec$wave
#' @param stLambda starting wavelength for fitting
#' @param endLambda ending wavelength for fitting
#' @param os1 tmp
#' @param os2 tmp
#' @param medWidth pixel  width of filtering kernal for MedSmooth 
#' @param smoothWidth pixel width of smoothing kernal for MedSmooth
#' @param runningWidth pixel width of smoothing kernal for RunningMeanSmooth
#' @param verbose TRUE/FALSE tell me what's going on
#' @author I. Baldry, L. Davies, L Drygala
#' @examples 
#' None applicable
#' @export
FitandFilter = function(spec, stLambda, endLambda, os1, os2, medWidth=51, smoothWidth= 121, runningWidth=21, verbose = TRUE){
  if(stLambda == -1){
    stLambda <- spec$lambda[1]
  } 
  if(endLambda == -1){
    endLambda <- tail(spec$lambda,1)
  }
  cFitRegion <- CriteriaBetween(spec$lambda, c(stLambda, endLambda))
  cFinite <- is.finite(spec$flux)
  useFit <- which(cFitRegion & cFinite)
  useNaN <- which(!cFinite)
  countNaN <- length(useNaN)
  # First fit ndegree polynomial and subtract fit. 
  # This uses an iterative rejection routine to remove outliers. 
  coeffs <- PolyReject(x = spec$lambda[useFit], y = spec$flux[useFit], verbose = verbose)
  
  specFit <- PolyCalc(coeffs, spec$lambda)
  spec$flux <- spec$flux - specFit
  spec$flux[useNaN] <- 0.0
  
  # Create smooth version of resulting spectrum - using median filtering --
  # and then a trapezium smoothing effectively - 
  # and then subtract smoothed spectrum. 
  fluxSmooth <- spec$flux
  #if numReject GE 1 then fluxSmooth[use_fit[reject_pts]] = 0.0
  fluxSmooth <- MedSmooth(fluxSmooth, medWidth, smoothWidth) #### old
  fluxSmooth <- RunningMeanSmooth(x = fluxSmooth, width = runningWidth) #### old
  spec$flux <- spec$flux - fluxSmooth
  
  # Set end points smoothly to zero between stLambda+os2 and stLambda+os1
  use <- which(spec$lambda < stLambda+os1)
  countStLambda <- length(use)
  spec$flux[use] <- 0.0
  use <- which(CriteriaBetween(spec$lambda, c(stLambda+os1,stLambda+os2)))
  count <- length(use)
  spec$flux[use] <- spec$flux[use] * CosineFilter(count,reverse = FALSE)
  
  # Set end points smoothly to zero between endLambda-os2 and endLambda-os1
  use <- which(spec$lambda > endLambda-os1)
  countEndLambda <- length(use)
  spec$flux[use] <- 0.0
  use <- which(CriteriaBetween(spec$lambda, c(endLambda-os2,endLambda-os1)))
  count <- length(use)
  spec$flux[use] <- spec$flux[use] * CosineFilter(count, reverse = TRUE)
  
  # Reset bad regions again to zero. 
  spec$flux[useNaN] = 0.0  
  return = list(spec, useFit)
}
