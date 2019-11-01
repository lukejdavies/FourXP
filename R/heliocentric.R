#' Calculate helocentric correction.
#'
#' @description Function compute correction term to add to velocities to 
#' convert to heliocentric reference frame. Either jd or tai must be supplied. 
#' 
#' @param ra Right ascension in degrees
#' @param dec Declination in degrees
#' @param epoch Epoch of observation for RA, DEC. Default is J2000.
#' @param jd Julian date
#' @param tai Number of seconds since Nov 17 1858
#' @param longitude Longitude of observatory
#' @param latitude Latitude of observatory
#' @param altitude altitude of observatory
#' @return Velocity correction term, in km/s, to add to measured radial velocity 
#' to convert it to the heliocentric frame.  
#' @author S. Burles
#' @author D. Schlegel
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' #Example for D10 field observed from the AAT on 18/12/2017.
#' vcorr<-Heliocentric(150.0, 1.55, epoch = 2000.0, tai = 58105.5, longitude =149.0661, latitude = -31.27704, altitude = 1164)
#' @export
Heliocentric = function(ra, dec, epoch = 2000.0, jd = FALSE, tai = FALSE, longitude = 360 - 105.820417, 
                        latitude = 32.780361, altitude = 2788){

  
  if (!jd && tai){
    jd <- 2400000.0 + tai / (24*3600)
  } else if (!jd){
    cat('\n *** Heliocentric ***Must specify either JD or TAI\n')
    return = 0
  }

  
  DRADEG <- 180 / pi
  
  #----------
  # Compute baryocentric velocity
  
  r <- baryvel(jd, epoch)
  dvelh <- r$dvelh
  dvelb <- r$dvelb

  # Project velocity toward star
  vbarycen <- dvelb[1]*cos(dec/DRADEG)*cos(ra/DRADEG) + 
    dvelb[2]*cos(dec/DRADEG)*sin(ra/DRADEG) + dvelb[3]*sin(dec/DRADEG) 
  
  #----------
  # Compute rotational velocity of observer on the Earth
  
  # LAT is the latitude in radians.
  latrad <- latitude / DRADEG
  
  # Reduction of geodetic latitude to geocentric latitude (radians).
  # DLAT is in arcseconds.
  
  dlat    <- -(11 * 60 + 32.743000) * sin(2 * latrad) + 
    1.163300 * sin(4 * latrad) -0.002600 * sin(6 * latrad)
  latrad  <- latrad + (dlat / 3600) / DRADEG
  
  # R is the radius vector from the Earth's center to the observer (meters).
  # VC is the corresponding circular velocity
  # (meters/sidereal day converted to km / sec).
  # (sidereal day = 23.934469591229 hours (1986))
  
  r <- 6378160 * (0.998327073 + 0.00167643800 * cos(2 * latrad) - 
                    0.00000351 * cos(4 * latrad) + 0.000000008 * cos(6 * latrad)) + altitude
  vc <- 2 * pi * (r / 1000)  / (23.934469591229 * 3600)
  
  # Compute the hour angle, HA, in degrees
  LST = ct2lst(longitude, 0, jd) # TODO 0 to make three params?

  LST <- 15 * LST # convert from hours to degrees
  HA <- LST - ra
  
  # Project the velocity onto the line of sight to the star.
  vrotate <- vc * cos(latrad) * cos(dec/DRADEG) * sin(HA/DRADEG)
  
  return = (vbarycen + vrotate)
  ##   return, vbarycen 
}                                                      