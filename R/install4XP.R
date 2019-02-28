#' Install/Download libraries for 4XP 
#'
#' @description Downlaod and unpack all things required for 
#' 
#' @param downloadModels TRUE/FALSE download models, if FALSE calibrators will just be unpacked
#' @param ModelDir where to unpack FourXP models
#' @return SpecDir (this will need to be used as a variable in FourXP_Sim() and makeSpec())
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' install4XP()
#' @export
install4XP<-function(downloadModels=TRUE, ModelDir='.'){
  

  for (i in 1:length(.libPaths())){
    tmp<-list.files(path=.libPaths()[i], pattern='*')
    if ('FourXP' %in% tmp ==TRUE){LibPath<-.libPaths()[i]}
  }
  
  system(paste('tar -xvf ', LibPath, '/FourXP/data/calibrators.tar --directory ', LibPath, '/FourXP/data/',sep='')) 

  
  if (downloadModels==TRUE){
    
    
    cat('** WARNING ** You need ~1Gb of free disk space to install 4XP model libraries. Do you wish to continue?', '\n')
    s<-readline('Continue (y/n):')
    if (s!='y'){
      cat('Stopping installation process', '\n')
      return(NULL)
    }
    
    
    cat('Installing FourXP models....', '\n')
    tmp<-list.files(path=ModelDir, pattern='*')
    if ('FourXPmodels.zip' %in% tmp ==TRUE){
      cat('Looks Like you already have the FourXPmodels.zip file in', ModelDir, 'do you wish to re-download it?', '\n')
      s<-readline('Re-download (y/n):')
      if (s!='y'){
        cat('Stopping installation process', '\n')
        return(NULL)
      }
    }
    if ('FourXPmodels' %in% tmp ==TRUE){
      cat('Looks Like you already have the FourXPmodels directory in', ModelDir, 'do you wish to re-download it?', '\n')
      s<-readline('Re-download (y/n):')
      if (s!='y'){
        cat('Stopping installation process', '\n')
        return(NULL)
      }
    }
    download.file('https://www.dropbox.com/s/hwnbabirn59796x/FourXPmodels.zip?raw=1',destfile=paste(ModelDir,'/FourXPmodels.zip',sep=''), method="auto")
    system(paste('unzip ', ModelDir,'/FourXPmodels.zip',sep=''))
    
    cat('\n', 'Download complete and data unpackaed')
    
    cat('\n\n', '4XP models are unpacked in ', paste(ModelDir,'/FourXPmodels/',sep=''), '- Use this location as specDir variable in FourXP_Sim() and makeSpec()', '\n\n')
    
  }
  
  return(paste(ModelDir,'/FourXPmodels/',sep=''))
  
  
}