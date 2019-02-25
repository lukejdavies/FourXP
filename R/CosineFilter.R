#' CosineFilter
#' @description Rising or falling cosine filter, i.e., a cosine bell taper.
#' @param num number of cosine filter
#' @param reverse TRUE/FALSE reverse filter
#' @author I. Baldry, L. Davies, L Drygala
#' @examples 
#' None applicable
#' @export
CosineFilter = function(num, reverse=FALSE){
  if(num<2){
    return = 0
  }
  phase <- (0:(num-1))  /(num-1) * pi
  if (!reverse){
    filter <- 0.5 * (1.0 - cos(phase)) 
  } else { 
    filter <- 0.5 * (1.0 + cos(phase))
  }
  return = filter
}
