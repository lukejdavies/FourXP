#' CriteriaBetween
#' @description Return 1 where data is within range
#' @param data 1D vector of data
#' @param data_range data range to return 1 within.
#' @author I. Baldry, L. Davies, L Drygala
#' @examples 
#' None applicable
#' @export
CriteriaBetween = function(data, data_range){
  criteria <- (data > data_range[1]) & (data <= data_range[2])
  return = criteria
}