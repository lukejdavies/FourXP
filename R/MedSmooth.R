#' MedSmooth
#' @description # Median filter a spectrum then smooth.
#' @param spec 1D vector to be filtered and smoothed
#' @param medWidth width of filtering kernal
#' @param smoothWidth width of smoothing kernal
#' @author I. Baldry, L. Davies, L Drygala
#' @examples 
#' None applicable
#' @export
MedSmooth = function(spec, medWidth = 51, smoothWidth = 121) {
  # Smooth with median filter. 
  outSpec <- MedianAdjust(spec,medWidth)
  
  # Further smooth using standard boxcar smoothing. 
  outSpec <- RunningMeanSmooth(x = outSpec, width = smoothWidth)
  
  return = outSpec
}
