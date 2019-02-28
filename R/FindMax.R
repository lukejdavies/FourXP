# Find maximum of yArray and fit quadratic using nearby points.
# Written by Ivan Baldry. 
# Translated by to R Leon Drygala
FindMax = function(xArray, yArray, n = 2){
  length <- length(xArray)
  xPts <- 1:length
  # find highest point of array, integer wise
  maxY <- max(yArray)
  xMaxPt = which(yArray == maxY)[[1]]
  #[[1]] ignores multiple maximums
  yAdjusted <- yArray - maxY
  # xAdjustment is a hack to make lm work properly on small x scales
  # TODO fix hack
  xAdjustment <- xArray[1]
  xAdjusted <- xArray - xAdjustment
  
  # fit quadratic to points +-n of xMaxPt
  usePts <- which( (xPts >= xMaxPt-n) & (xPts <= xMaxPt+n) )
  result <- lm(yAdjusted[usePts] ~ poly(xAdjusted[usePts], 2, raw = TRUE))
  if (result$rank == 3) {
    # Solve for dY/dX <- 0:  2*result$coefficients[[3]]*X + result$coefficients[[2]] = 0
    xCenter = -0.5 * result$coefficients[[2]] / result$coefficients[[3]]
  } else {
    cat('\nFindMax: quadratic fit failed\n')
    xCenter = xArray[xMaxPt]
  }
  return = xCenter + xAdjustment
}
