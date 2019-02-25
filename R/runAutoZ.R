#' Run AutoZ over a supplied list of spectra
#'
#' @description This is the highlevel main TAZ function for running AutoZ over a number of 
#' spectra, or all spectra in the 'data/reduced/stackedSpec/' directory. Function will run AutoZ
#' and make plots in the 'data/reduced/stackedSpec/AutoZplots' directory.
#'  
#' @param specs vector list of file paths to run AutoZ over. Must either be full directory path, or can
#' set to 'all' to run over all spectra in the 'data/reduced/stackedSpec/' directory.  
#' @param logName log filename to write progress to
#' @param verbose tell me whats going on: 0=nothing, 1=somethings, 2=everything
#' @param makePlots TRUE/FLASE make plots of spectra using plotSpec.R
#' @param cores number of cores to use
#' @examples 
#' runAutoZ(specs='all', logName='tempLog.txt', verbose=1, makePlots=T, cores=4)
#' 
#' runAutoZ(specs=c('data/reduced/stackedSpec/G007372.RData'), logName='tempLog.txt', verbose=1, makePlots=T, cores=4)
#' @export
runAutoZ<-function(specs=specs, logName=logName, verbose=verbose, makePlots=T, cores=cores, highZ=T){


    if (specs[1]=='all'){
        specs<-c()
        
        listSpec<-list.files(path='data/reduced/stackedSpec/', pattern='*.Rdata')
        specs<-paste('data/reduced/stackedSpec/', listSpec, sep='')
        
    }
    
    
    registerDoParallel(cores=cores)

    a = foreach(i=1:length(specs)) %dopar%  {
    #for (i in 1:length(specs)){

        if (verbose>1){cat('    - Running AutoZ for spectrum: ',specs[i],  '\n')}
        write(paste('    - Running AutoZ for spectrum: ',specs[i],sep=''), file=logName, append=T)

        load(specs[i])
        
        spec2<-spec
        spec2$flux<-spec$fluxSub
        spec2$error<-spec$sn
        spec2$longitude =149.0661
        spec2$latitude = -31.27704
        spec2$altitude = 1164
        spec2$error<-sqrt(spec2$error)

        
        autoz_out<-AutozSingleSpec(spec2,tempFile = 'data/calibrators/AutoZTemp/filtered-templates.fits',doHelio=F, verbose=F,highZ=highZ, templateNumbers = c(2:14,16:22,40:47))
        
        spec$z<-autoz_out$results[1]
        spec$prob<-autoz_out$prob
        spec$cc<-autoz_out$results[2]
        spec$z2<-autoz_out$results[4]
        spec$cc2<-autoz_out$results[5]
        spec$Temp<-autoz_out$results[3]

        if (verbose>1){cat('        - AutoZ found redshift of ',spec$z,  '\n')}
        write(paste('        - AutoZ found redshift of ',spec$z, sep=''), file=logName, append=T)
        if (verbose>1){cat('        ...with probability of ',spec$prob,  '\n')}
        write(paste('        ...with probability of : ',spec$prob,sep=''), file=logName, append=T)
 
        save(spec,file=specs[i])

          if (makePlots==T){
  
              if (verbose>1){cat('        - Plotting AutoZ outputs as: ',paste('data/reduced/stackedSpec/AutoZplots/', spec$ID,'.pdf',sep=''),  '\n')}
              write(paste('        - Plotting AutoZ outputs as: data/reduced/stackedSpec/AutoZplots/', spec$ID,'.pdf', sep=''), file=logName, append=T)
              
              pdf(paste('data/reduced/stackedSpec/AutoZplots/', spec$ID,'.pdf',sep=''), width=18, height=18)
              
              plotSpec(spec)
              
              dev.off()
  
              }
          
    }

}
