#' Compute cross corelations in FourXP_Z
#'
#' @description Function to compute cross-correlation strngth between a given spectrum and
#' a set of templates. Function is internal to FourXP_Z.
#' 
#' @param spec input spectrum containing spec$lambda and spec$flux
#' @param gap gap in shifts probed in cross correlation
#' @param tempDatar structure of template data 
#' @param heloval helocentric correction value
#' @param plan fast-forier tranform plan
#' @param z_prior redshift priors, c(lo,hi)
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' None applicable as internal function.... 
#' @export
DoCrossCorr = function(spec, gap, tempData, helioVel, plan, z_prior, highZ=FALSE){
  
  corrSize <- 17000
  if (highZ==T){corrSize <- 24999}
  
  # prepare output data
  dataout <- vector('list',length(tempData$templateNumbers))
  count = 1
  #plan = planFFT(length(spec$flux), effort = 2)

  for( tempNum in tempData$templateNumbers ){
    # Find index of certain template number                                         
    tempID <- which(tempData$templateNumbers == tempNum)
#takes about 0.224 to next section
    crossCorrRaw <- CrossCorr(template = tempData$specData[[tempID]]$flux, spec = spec$flux, 
                              plan = plan, corrSize = corrSize)
#takes about 0.046 to next section
    # calculate redshifts
    shifts <- 10^((( 0:(2*corrSize) ) - corrSize) * gap) *
      (1 + tempData$specData[[tempID]]$redshift) * (1 + helioVel/2.998E5)   - 1.0
    
    
#takes about 0.002 to next section
    r <- GetRange(tempNum, z_prior)
    rmsZRange     <- r[[1]]
    allowedZRange <- r[[2]]
#takes about 0.048 to next section
    # set criteria
    criteriaSearch  <- CriteriaBetween(shifts, allowedZRange)
    criteriaRMS     <- CriteriaBetween(shifts, rmsZRange)
#takes about 0.214 to next section
    criteriaPospeak <- CriteriaPeak(crossCorrRaw)  # positive peaks
    criteriaNegpeak <- CriteriaPeak(-crossCorrRaw) # negative peaks
# takes about 0.028  to next section
    criteriaTP      <- criteriaPospeak | criteriaNegpeak
# takes about 0.049  to next section
    # Subtract trimmed mean excluding top and bottom 4% of points.
    # This brings more symmetry to positive and negative peaks.
    useRMS    <- which(criteriaRMS)
    countRMS  <- length(useRMS)
    crossCorr <- crossCorrRaw - MeanReject(crossCorrRaw[useRMS], countRMS/25)
# takes about 0.034 to next section
    # normalisation using turning points - divide by root mean square.
    useNorm <- which(criteriaRMS & criteriaTP)
    numTurningpoints <- length(useNorm)
    crossCorr <- crossCorr / sqrt( mean((crossCorr[useNorm])^2 ) )
# this all takes about 0.081 to next section
    # TODO this is commented out in Ivans code, not sure to include
    # normalization using values of positive peaks and 
    # negative values of negative peaks.
    usePos            <- which(criteriaRMS & criteriaPospeak)
    countPos          <- length(usePos)
    useNeg            <- which(criteriaRMS & criteriaNegpeak)
    countNeg          <- length(useNeg)
    numTurningpoints  <- countPos + countNeg
    testVals          <- c(crossCorr[usePos], -crossCorr[useNeg])
    trimmedMean       <- MeanReject(testVals, numTurningpoints/25)
    sdEstimate        <- rms(testVals - trimmedMean)
    crossCorr         <- (crossCorr - trimmedMean) / sdEstimate

# this take 0.041 seconds to next section
    # assign information to structure for output function
    maskInfo <- 1*criteriaSearch + 2*criteriaRMS + 4*criteriaPospeak + 
                8*criteriaNegpeak + 16*criteriaTP
    
    dataout[[count]] <- list("templateNumber" = tempNum,"shifts"=shifts, "crossCorrRaw" = crossCorrRaw,
                            "crossCorr"=crossCorr,"maskInfo" = maskInfo,"numTurningpoints" = numTurningpoints)
    count <- count + 1
  }
#time <- proc.time() 
#cat("\nCROSS_CORR totalCorrTime is:", totalCorrTime)

  return = dataout
}
