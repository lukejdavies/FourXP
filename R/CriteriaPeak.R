#' CriteriaPeak
#'
#' @description Determine where local peak values are in data. 
#' @param values 1D data values
#' @author I. baldry, L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' None applicable as internal function.... 
#' @export
CriteriaPeak = function(values){
  num <- length(values)
  out <- rep(FALSE, num)
  out[2:(num-3)] =  (values[2:(num-3)] >= values[1:(num-4)]) & 
    (values[2:(num-3)] > values[3:(num-2)])
  return = out
}
