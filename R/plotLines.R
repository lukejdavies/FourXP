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
#' @param AGNCol colour to plot strong AGN lines
#' @param lty line type for plotting
#' @param lwd line thinckness for plotting
#' @param cex label text size
#' @param main TRUE/FALSE only plot common lines
#' @param Air TRUE/FALSE plot line in air wavelengths (FALSE=Vacuum wavelengths)
#' @author L. Davies <luke.j.davies@uwa.edu>

#' @examples
#' load(paste(.libPaths(),'/FourXP/data/ExampleSpec.Rdata',sep=''))
#' plot(spec$wave, hanning.smooth(spec$flux, degree=9), type='l', xlab='Wavelength, ang', ylab='Counts') 
#' plotLines(z=spec$z, abPos=100, labOff=-50, EmCol='blue', AbsCol='darkgreen', lty=2, lwd=2, cex=0.5)
#' @export
plotLines<-function(z=0, xunit='ang', labPos=100, labOff=0, EmCol='blue', AbsCol='darkgreen', AGNCol='skyblue2', lty=2, lwd=1, cex=0.3,main=F, Air=TRUE, ...){


    lines<-load.lines()
    
    if (Air==T){lines$wave_ang<-AirFromVacuum(lines$wave_ang)}
    
    if (main==T){lines<-lines[which(lines$names %in% c('Ly-alpha', 'O II', '', 'H-beta', 'O III[4932]', '','O III[5008]', 'N II', 'H-alpha', 'N II', 'S II', 'S II', 'CaK&H', '', 'G-band', 'Mg', 'Na')==T),]}
    
    if (xunit=='ang') {line_x <- as.numeric(lines$wave_ang)*(1+z)}
    if (xunit=='hz') {line_x <- as.numeric((2.988e8)/(lines$wave_ang/10^10))/(1+z)}
    if (xunit=='micron') {line_x <- as.numeric(lines$wave_ang/10^4)*(1+z)}
    if (xunit=='m') {line_x <- as.numeric(lines$wave_ang/10^10)*(1+z)}
    if (xunit=='nm') {line_x <- as.numeric(lines$wave_ang/10)*(1+z)}

    
    
    for (i in 1:length(lines$names)) {
        if (lines$stellar[i]==F & lines$agnStrong[i]==F) {
            
            #abline(v=line_x[i], col=EmCol, lty=lty, lwd=lwd)
          
            lines(c(line_x[i],line_x[i]), c(par("usr")[3],labPos-(par("usr")[4]-par("usr")[3])/40), col=EmCol, lty=lty, lwd=lwd)
            
            if (lines$names[i] %in% c('Ly-alpha','H-delta','H-gamma', 'H-beta', 'H-alpha')==F){text(line_x[i]+labOff, labPos, lines$names[i], col=EmCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='Ly-alpha'){text(line_x[i]+labOff, labPos, expression('Ly'~alpha), col=EmCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='H-delta'){text(line_x[i]+labOff, labPos, expression('H'~delta), col=EmCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
          if (lines$names[i]=='H-gamma'){text(line_x[i]+labOff, labPos, expression('H'~gamma), col=EmCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
          if (lines$names[i]=='H-beta'){text(line_x[i]+labOff, labPos, expression('H'~beta), col=EmCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
          if (lines$names[i]=='H-alpha'){text(line_x[i]+labOff, labPos, expression('H'~alpha), col=EmCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
        }
        if (lines$stellar[i]==T & lines$agnStrong[i]==F) {
          lines(c(line_x[i],line_x[i]), c(par("usr")[3],labPos-(par("usr")[4]-par("usr")[3])/40), col=AbsCol, lty=lty, lwd=lwd)
            text(line_x[i]+labOff, labPos, lines$names[i], col=AbsCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))
        }
          if (lines$agnStrong[i]==T) {
            lines(c(line_x[i],line_x[i]), c(par("usr")[3],labPos-(par("usr")[4]-par("usr")[3])/40), col=AGNCol, lty=lty, lwd=lwd)
            if (lines$names[i] %in% c('Ly-alpha','H-delta','H-gamma', 'H-beta', 'H-alpha')==F){text(line_x[i]+labOff, labPos, lines$names[i], col=AGNCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='Ly-alpha'){text(line_x[i]+labOff, labPos, expression('Ly'~alpha), col=AGNCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='H-delta'){text(line_x[i]+labOff, labPos, expression('H'~delta), col=AGNCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='H-gamma'){text(line_x[i]+labOff, labPos, expression('H'~gamma), col=AGNCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='H-beta'){text(line_x[i]+labOff, labPos, expression('H'~beta), col=AGNCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
            if (lines$names[i]=='H-alpha'){text(line_x[i]+labOff, labPos, expression('H'~alpha), col=AGNCol,cex=cex, srt=270, pos=2, offset=c(0.0,0.5))}
          }
    }
    
   
    
          


}

