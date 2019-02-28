# Calculate figure of merit and and certainty of redshift
# Translated from I. Baldry AutoZ code
# Written by Leon Drygala
CalculateCertainty = function(spec){
  # Calculate root mean square to mean absolute deviation ratio. 
  spec$rmsMadRatio   <- spec$rmsNorm / spec$meanadNorm
  
  # Calculate ratio of first peak compared to RMS of 2nd/3rd/4th. 
  spec$ccSigma1to234 <- spec$ccSigma[1] / rms(spec$ccSigma[2:4])
  
  # The following calibrations were determined for the GAMA survey. 
  # It is unclear how they will perform with different surveys. 
  
  # Define FOM, lower of cc_sigma and a converted cc_sigma1to234.
  spec$ccFOM <- PolyCalc(c(0.4, 2.8), spec$ccSigma1to234)[1]
  if(spec$ccFOM > spec$ccSigma[1]) spec$ccFOM <- spec$ccSigma[1]
  
  # Adjustment to FOM for large rms_mad_ratio. 
  rmrAdjustment <- (1.5 * (spec$rmsMadRatio - 1.8)) 
  if(rmrAdjustment < 0) rmrAdjustment <- 0.0
  if(rmrAdjustment > 2.1) rmrAdjustment <- 2.1
  if(spec$ccFOM - rmrAdjustment < 2.6 ){
    spec$ccFOM <- 2.6
  } else {
    spec$ccFOM <- (spec$ccFOM - rmrAdjustment)
  }
  
  # Set redshift confidence from CC_FOM. 
  x = (spec$ccFOM - 3.70) / 0.70
  spec$prob <- (tanh(x) + 1) / 2
  
  return = spec
}