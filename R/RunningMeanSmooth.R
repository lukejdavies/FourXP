#' RunningMeanSmooth
#' @description Translation of IDL's SMOOTH function into R Uses IDL's '/EDGE_TRUNCATE' option
#' @param x 1D vector to be smoothed
#' @param width width of smoothing kernal
#' @author L Drygala
#' @examples 
#' None applicable
#' @export
RunningMeanSmooth = function(x,width){
  length <- length(x)
  #filter doesn't seem to be able to deal with edge values
  outX <- filter(x, rep(1/width,width), method = "convolution", sides = 2)
  
  # Truncate left edge
  for(i in 1:(width%/%2)){
    index <- (i-width%/%2):(i+width%/%2)
    for(j in 1:width){
      if(index[j]<1) index[j] <- 1
    }
    outX[i] <- mean(x[index])
  }
  
  # Truncate right edge
  for(i in (length - width%/%2 + 1):length){
    index <- (i-width%/%2):(i+width%/%2)
    for(j in 1:width){
      if(index[j]>length) index[j] <- length
    }
    outX[i] <- mean(x[index])
  }
  
  return = outX
}
