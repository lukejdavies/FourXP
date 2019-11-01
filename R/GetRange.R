#' GetRange
#'
#' @description Set allowable redshfit ranges of AutoZ Templates
#' 
#' @param tNum input spectrum containing spec$lambda and spec$flux
#' @param z_prior user-define photo-z prior
#' @author I. baldry, L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' None applicable as internal function.... 
#' @export
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
  } else if (tNum == 29 | tNum == 32){
    # QSO templates - not working reliably with GAMA - needs highz set
    rmsZRange <- c(-0.1,5)
    allowedZRange <- c(0.8, 5.500)
  } else if (tNum == 30 | tNum == 31){
    # QSO templates - not working reliably with GAMA - needs highz set
    rmsZRange <- c(-0.1,5)
    allowedZRange <- c(1.5, 5.500)
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
