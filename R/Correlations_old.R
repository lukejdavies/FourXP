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
#' @author I. baldry, L. Davies <luke.j.davies@uwa.edu>
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

# Cross correlate one spectrum with a template
# "template" and "spec" must be the same length (even number). 
# Return a reduced size of array centred around zero shift. 
# Written ages ago in IDL by Ivan Baldry. 
# Translated to R by Leon Drygala
CrossCorr = function(template, spec, plan, corrSize = 14000){
  
  tempLength = length(template)

  # take conjugate of fourier transform of rebinned template spectra
  #fftTemplate <- FFT(template, plan = plan) / tempLength#TODO consider other (faster) transforms
  fftTemplate <- fft(template)/ tempLength
  fftTemplate <- Conj(fftTemplate)
  
  #take fourier transform of spectra
  #fftSpec <- FFT(spec, plan = plan) / tempLength
  fftSpec <- fft(spec) / tempLength
  
  # multiply by conj. of fourier transform of t_plate spectrum
  fftSpec <- fftSpec * fftTemplate
  #invfft <- IFFT(fftSpec, plan = plan, scale = FALSE)
  invfft <- fft(fftSpec,inverse=TRUE)
  
  #take real part of inverse FFT for cross correlation. 
  crC <- Re(invfft) 
  length <- length(crC)
  halfLength <- length/2
  
  crC <- c(crC[(halfLength+1):(length)], crC[1:(halfLength)])
  
  #create output array: length = 2*corr_size+1
  crossCorr <- crC[(halfLength-corrSize+1):(halfLength+corrSize+1)]
  return = crossCorr
}



# Determine where local peak values are in data. 
# Written by Ivan Baldry. 
# Translated by Leon Drygala
CriteriaPeak = function(values){
  num <- length(values)
  out <- rep(FALSE, num)
  out[2:(num-3)] =  (values[2:(num-3)] >= values[1:(num-4)]) & 
                    (values[2:(num-3)] > values[3:(num-2)])
  return = out
}

# Determine where local peak values are in data. 
# Attempt to increase performance over original CriteriaPeak
# Returns 1 for neg peak 2 for pos peak 0 for no peak
# Writen by Leon Drygala
CriteriaPeakFast = function(values){
  num <- length(values)
  out <- vector(mode = 'logical',length = num)
  
  out[2:(num-1)] = ( (values[2:(num-1)] >= values[1:(num-2)]) == 
    (values[2:(num-1)] > values[3:(num)]) ) 
  
  out[out] <- out[out] + (values[out] > values[c(out[-1],F)])
  
  return = out
}


# Take the mean of a set of values after rejecting numReject lowest 
# and numReject highest values.
# This is called the 'trimmed mean' or 'truncated mean'.
# Written by Ivan Baldry. 
# Translated by Leon Drygala
MeanReject = function(data, numReject){
  if (numReject){
    # sort data
    dataValues <- sort(data)
    num <- length(dataValues)
    
    stPoint <- numReject + 1
    endPoint <- num - numReject
    
    result <- mean(dataValues[stPoint:endPoint])
  } else {
    result <- mean(data)
  }
  return = result
}

# Take the mean of a set of values after rejecting numReject lowest 
# and numReject highest values.
# This is called the 'trimmed mean' or 'truncated mean'.
# Rewritten to try and improve speed and avoid sorting
# Written by Leon Drygala
MeanRejectFast = function(data, numReject){
  tooLarge <- PriorityQueue(decreasing = F)
  tooSmall <- PriorityQueue(decreasing = T)
  for(i in 1:numReject){
    tooLarge[['insert']](data[i], i)
    tooSmall[['insert']](data[i], i)
  }
  
  for(i in (numReject+1):length(data)){
    if(data[i] > tooLarge[['peak']]()[1]){
      tooLarge[['insert']](data[i], i)
      tooLarge[['pop']]()
    }
    if(data[i] < tooSmall[['peak']]()[1]){
      tooSmall[['insert']](data[i], i)
      tooSmall[['pop']]()
    }
  }
  rejectIndex = c(tooSmall[['dump']]()[[2]], tooLarge[['dump']]()[[2]])
  newData = data[-rejectIndex]
  result = mean(newData)
  return = newData
}

PriorityQueue <- function(decreasing = T) {
  keys <<- values <<- NULL
  insert <- function(key, value) {
    temp <- c(keys, key)
    ord <- order(temp, decreasing = decreasing)
    keys <<- temp[ord]
    values <<- c(values, c(value))[ord]
  }
  peak <- function() {
    head <- c(keys[[1]], values[[1]])
    return(head)
  }
  pop <- function() {
    head <- c(keys[[1]], values[[1]])
    values <<- values[-1]
    keys <<- keys[-1]
    return(head)
  }
  
  queueLength <- function() {
    return(length(values))
  }
  dump <- function(){
    return(list(keys,values))
  }
  empty <- function() length(keys) == 0
  list(insert = insert, pop = pop, empty = empty, dump = dump, 
       peak = peak, queueLength = queueLength)
}


GetRange = function(tNum, z_prior){
  #TODO what should default be?? templates 1->10, 16, 18, 20, 21 not covered
  rmsZRange <- c(-0.1,0.5)
  allowedZRange <- c(-0.002, 0.002)
  if (sum(tNum == c(11,12,13,14,15,17,19,22))){
    # late-type stellar templates - power is at red end
    rmsZRange <- c(-0.2,0.4)
    allowedZRange <- c(-0.002, 0.002)
  } else if (tNum<=22){
    # remaining stellar templates
    rmsZRange <- c(-0.1,0.5)
    allowedZRange <- c(-0.002, 0.002)
  } else if ( tNum >= 23 && tNum <= 28){
    # original galaxy templates
    rmsZRange <- c(-0.1,0.8)
    rmsZRange <- c(-0.1,1.5)
    allowedZRange <- c(-0.005, 1.500)
  } else if (tNum >= 29 && tNum <= 32){
    # QSO templates - not working reliably with GAMA - needs highz set
    rmsZRange <- c(-0.1,5)
    allowedZRange <- c(0.8, 5.500)
  } else if (tNum >= 33 && tNum <= 49){
    # other galaxy templates
    rmsZRange <- c(-0.1,0.9)
    rmsZRange <- c(-0.1,1.5)
    allowedZRange <- c(-0.005, 1.500)
  } else if (tNum >= 50 && tNum < 60){
    rmsZRange <- c(-0.1,2.0)
    allowedZRange <- c(-0.005, 2.000)
  }  else if (tNum >=60 && tNum <= 72 ){
    rmsZRange <- c(-0.1,5.)
    allowedZRange <- c(-0.005, 2.000)
    
  } else if (tNum >73 && tNum <= 80 ){
    rmsZRange <- c(-0.1,10.0)
    allowedZRange <- c(-0.005, 2.000)
    
  } else if (tNum > 80){
    rmsZRange <- c(-0.1,2.0)
    allowedZRange <- c(-0.005, 2.000) #TODO this should be = z_prior, look up in do_crossvorr.pro in IDL
  }
  
  if (tNum ==64) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==65) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==66) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==67) {allowedZRange <- c(2.0, 6.5)}  
  
  if (tNum ==76) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==77) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==78) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==79) {allowedZRange <- c(2.0, 6.5)}
  if (tNum ==80) {allowedZRange <- c(2.0, 6.5)}
  
  if (allowedZRange[1]<z_prior[1]){allowedZRange[1]<-z_prior[1]}
  if (allowedZRange[2]>z_prior[2]){allowedZRange[2]<-z_prior[2]}
  
  return = list(rmsZRange, allowedZRange)
}

rms = function(x){
  return = sqrt(sum(x^2)/length(x))
}

