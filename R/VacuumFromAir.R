#' Internal AutoZ routines for converting vacuum wavelengths to air wavelengths 
#'
#' @description Function interpolate air wavelengths to vacuum wavelengths
#' using a lookup table. Equation comes from the SDSS web page:
#' http://www.sdss.org/dr7/products/spectra/vacwavelength.html. 
#'  Reference to Morton (1991, ApJS, 77, 119). 
#' 
#' @param air Wavelength in air
#' @return Wavelength in vacuum
#' @author I. Baldry
#' @examples 
#' Vac<-VacuumFromAir(3727.456)
#' @export
VacuumFromAir = function(air){
  
  vacLookup <- 0:7999 + 3500
  convertLookup <- (1.0 + 2.735182E-4 + 131.4182 / vacLookup^2 + 2.76249E8 / vacLookup^4)
  airLookup <- vacLookup / convertLookup
  convert <- approx(x = airLookup, y = convertLookup, xout = air, method = "linear",rule=2)$y
  vac <- air * convert
  
  return = vac
}

#' Internal AutoZ routines for converting air wavelengths to vacuum wavelengths 
#'
#' @description Function interpolate vacuum wavelengths to air wavelengths
#' using a lookup table. Equation comes from the SDSS web page:
#' http://www.sdss.org/dr7/products/spectra/vacwavelength.html. 
#'  Reference to Morton (1991, ApJS, 77, 119). 
#' 
#' @param Vac Wavelength in vacuum
#' @return Wavelength in air
#' @author I. Baldry
#' @examples 
#' Air<-AirFromVacuum(3728.516)
#' @export
AirFromVacuum = function(vac){
  air <- vac / (1.0 + 2.735182E-4 + 131.4182 / vac^2 + 2.76249E8 / vac^4)
  return = air
}