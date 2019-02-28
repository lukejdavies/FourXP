# Find highest peaks across the cross-correlation functions. 
# Note that in this routine peaks.template and templateList refers
# to the index within the ccinfo structure and NOT the template number that
# is finally assigned. 
# Written by Ivan Baldry.
# Translated to R by Leon Drygalaa
FindHighestPeaks = function(ccinfo, num = 4){
  templateList    <- c(-1)
  templateIDList  <- c(-1)
  indexesList     <- c(-1)
  redshiftList    <- c(-1)
  crosscorrList   <- c(0)
  
  # first extract all peaks in normalised cross correlation value
  # across all templates
  for (i in 1:length(ccinfo)) {
    # cross correlation values of peaks in allowed range (not a peak = 0)
    #testvals <- ccinfo[[i]]$crossCorr * ((ccinfo[[i]]$maskInfo & 5) == 5)#TODO what on earth does this do?
    testvals <- ccinfo[[i]]$crossCorr * (ccinfo[[i]]$maskInfo == 23)#TODO I think this is what it's meant to do...
    # save only certain peaks
    use <- which(testvals >= 1)
    count <- length(use)
    if (count >= 1) {
      templateList    <- append(templateList, rep(0,count) + ccinfo[[i]]$templateNumber)
      templateIDList  <- append(templateIDList, rep(0,count) + i)
      indexesList     <- append(indexesList, use)
      redshiftList    <- append(redshiftList, ccinfo[[i]]$shifts[use])
      crosscorrList   <- append(crosscorrList, testvals[use])
    }
  }
  
  peaks <- vector('list',num)
  for (i in 1:num) {
    # find highest peak remaining
    maxcorr <- max(crosscorrList)
    pt = which(crosscorrList == maxcorr)[[1]]
    #[[1]] ignores any multiple maximums
    peaks[[i]] <- list("template" = templateList[pt], "redshift" = redshiftList[pt], 
                       "crossCorr"=maxcorr, "shiftIndex"=indexesList[pt], "templateID" = templateIDList[pt])
    
    # exclude +- 600 km/s in next search
    zrange  <- (1 + peaks[[i]]$redshift) * (1 + c(-0.002,0.002)) - 1.
    criteriaZRange <- (redshiftList > zrange[1]) & (redshiftList <= zrange[2])
    use <- which(criteriaZRange)
    crosscorrList[use] <- -1.
  }
  return = peaks
}
