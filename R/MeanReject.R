#' MeanReject
#'
#' @description Take the mean of a set of values after rejecting numReject lowest 
# and numReject highest values.This is called the 'trimmed mean' or 'truncated mean'.
#' @param data 1D data values
#' @param numReject number of value to reject
#' @author I. baldry, L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' None applicable as internal function.... 
#' @export
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
