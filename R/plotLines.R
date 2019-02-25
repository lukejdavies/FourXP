#' Overplot the position of common spectral features on the currently plotted spectrum 
#'
#' @description Function to plot the position of common spectral absorption and emission line features
#' over the current plotting window. 

#' @param z redshift at which to overplot lines
#' @param xunit units of the xaxis. Options are ang,hz,micron,m,nm
#' @param labPos y axis value of where to print line labels
#' @param labOff xaxis offset to apply to labels
#' @param EmCol colour to plot emission lines
#' @param AbsCol colour to plot absorption lines
#' @param lty line type for plotting
#' @param lwd line thinckness for plotting
#' @param cex label text size
#' @author L. Davies <luke.j.davies@uwa.edu>

#' @examples
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plot(spec$wave, hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts') 
#' plotLines(z=spec$z, abPos=100, labOff=-50, EmCol='blue', AbsCol='darkgreen', lty=2, lwd=2, cex=0.5)
#' @export
plotLines<-function(z=0, xunit='ang', labPos=100, labOff=-50, EmCol='blue', AbsCol='darkgreen', lty=2, lwd=1, cex=0.3,...){


    lines<-load.lines()
    
    if (xunit=='ang') {line_x <- as.numeric(lines$wave_ang)*(1+z)}
    if (xunit=='hz') {line_x <- as.numeric(lines$freq_hz)/(1+z)}
    if (xunit=='micron') {line_x <- as.numeric(lines$wave_micron)*(1+z)}
    if (xunit=='m') {line_x <- as.numeric(lines$wave_m)*(1+z)}
    if (xunit=='nm') {line_x <- as.numeric(lines$wave_nm)*(1+z)}

    
    
    for (i in 1:length(lines$names)) {
        if (lines$stellar[i]==F) {
            
            abline(v=line_x[i], col=EmCol, lty=lty, lwd=lwd)
            text(line_x[i]+labOff, labPos, lines$names[i], col=EmCol,cex=cex, srt=270)
        }
        if (lines$stellar[i]==T) {
            abline(v=line_x[i], col=AbsCol, lty=lty, lwd=lwd)
            text(line_x[i]+labOff, labPos, lines$names[i], col=AbsCol,cex=cex, srt=270)
        }
    }
    
   
    
          


}

