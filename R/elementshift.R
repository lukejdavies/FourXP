#' Apply element shift to vector in Hanning Smoothing for plotting
#'
#' @description Apply element shift to vector in Hanning smoothing function.
#' 
#' @param x vector
#' @param dx shift required

#' @examples 
#' None applicable as internal function.... 
#' @export
elementshift <-
function(x, dx) {
  if (dx<0) { dx<-length(x)-(abs(dx)%%length(x)) }
  return=(c(x[-(dx%%length(x)):0],x[1:(dx%%length(x))]))
}

