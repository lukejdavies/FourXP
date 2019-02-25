#' Reads in filter responses for a specified filter.
#'
#' @description Function reads in a specified filter tansmission curve from the FourXP data directory. Note that for this to
#' work you must have unpacked the FourXP data files using:
#' 
#' @description > LibPaths<-.libPaths()
#' @description > system(paste('tar -xvf ', LibPaths, '/FourXP/data/calibrators.tar --directory ', LibPaths, '/FourXP/data/',sep='')) 
#' 
#' @author L. Davies <luke.j.davies@uwa.edu>
#' @param filter string of the filter to load. Options are:
#' GALEX_FUV, GALEX_NUV, sdss_u,sdss_g,sdss_r,sdss_i,sdss_z, des_u,des_g,
#' des_r,des_i,des_z,des_Y,VISTA_Z,VISTA_Y,VISTA_J,VISTA_H,VISTA_K, WISE_W1,
#' WISE_W2,WISE_W3,WISE_W4,Herschel_100,Herschel_160,Herschel_250,Herschel_350,Herschel_500
#' 
#' @return the filter response in the form c(rowNum, wavelength(ang), transmission)
#' @examples 
#' NUV_filt<-getfilt('GALEX_NUV')
#' u_filt<-getfilt('sdss_u')
#' r_filt<-getfilt('sdss_r')
#' Z_filt<-getfilt('VISTA_Z')
#' J_filt<-getfilt('VISTA_J')
#' K_filt<-getfilt('VISTA_K')
#' 
#' plot(r_filt[,2], r_filt[,3], type='l', col='green', xlab='Wavelength, Ang', ylab='Transmission', xlim=c(1000,25000), ylim=c(0,1))
#' lines(NUV_filt[,2], NUV_filt[,3], col='purple')
#' lines(u_filt[,2], u_filt[,3], col='cyan')
#' lines(Z_filt[,2], Z_filt[,3], col='orange')
#' lines(J_filt[,2], J_filt[,3], col='red')
#' lines(K_filt[,2], K_filt[,3], col='brown')
#' @export
getfilt=function(filter){
  out=NA
  
  LibPaths<-.libPaths()
  
  FourXP_PATH<-NA
  
  for (jj in 1:length(LibPaths)){
    packagesList<-list.files(path=LibPaths[jj], pattern='*')
    if (length(which(packagesList=='FourXP'))>0){FourXP_PATH<-LibPaths[jj]}
  }
  
  
  if(filter=='GALEX_FUV'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_FUV_GALEX.tab', sep=''))     
  }
  if(filter=='GALEX_NUV'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_NUV_GALEX.tab', sep=''))     
  }
  if(filter=='sdss_u'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_u_SDSS.tab', sep=''))     
  }
 if(filter=='sdss_g'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_g_SDSS.tab', sep=''))     
 }
   if(filter=='sdss_r'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_r_SDSS.tab', sep=''))     
   }
   if(filter=='sdss_i'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_i_SDSS.tab', sep=''))     
   }
   if(filter=='sdss_z'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_z_SDSS.tab', sep=''))     
   }

  
  if(filter=='des_u'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_u_DES.tab', sep=''))     
  }
 if(filter=='des_g'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_g_DES.tab', sep=''))     
 }
   if(filter=='des_r'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_r_DES.tab', sep=''))     
   }
   if(filter=='des_i'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_i_DES.tab', sep=''))     
   }
   if(filter=='des_z'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_z_DES.tab', sep=''))     
   }

     if(filter=='des_Y'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_Y_DES.tab', sep=''))     
     }

     if(filter=='VST_u'){
       out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_u_VST.tab', sep=''), header=T)
       out<-cbind(seq(1,length(out[,1]),1), out[,1], out[,2])
  }
 if(filter=='VST_g'){
     out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_g_VST.tab', sep=''), header=T)
     out<-cbind(seq(1,length(out[,1]),1), out[,1], out[,2])
 }
   if(filter=='VST_r'){
       out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_r_VST.tab', sep=''), header=T)
       out<-cbind(seq(1,length(out[,1]),1), out[,1], out[,2])
   }
   if(filter=='VST_i'){
       out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_i_VST.tab', sep=''), header=T)
       out<-cbind(seq(1,length(out[,1]),1), out[,1], out[,2])
   }
   if(filter=='VST_z'){
       out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_z_VST.tab', sep=''), header=T)
       out<-cbind(seq(1,length(out[,1]),1), out[,1], out[,2])
   }


  
   if(filter=='VISTA_Z'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_Z_VISTA.tab', sep='')) 
      out[,3]<-out[,3]/100
   }
   if(filter=='VISTA_Y'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_Y_VISTA.tab', sep=''))
      out[,3]<-out[,3]/100
  }

   if(filter=='VISTA_J'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_J_VISTA.tab', sep='')) 
      out[,3]<-out[,3]/100
  }

   if(filter=='VISTA_H'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_H_VISTA.tab', sep='')) 
      out[,3]<-out[,3]/100
  }

   if(filter=='VISTA_K' | filter=='VISTA_Ks'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_K_VISTA.tab', sep=''))
      out[,3]<-out[,3]/100
   }
   if(filter=='WISE_W1'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_W1_WISE.tab', sep=''))     
   }
   if(filter=='WISE_W2'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_W2_WISE.tab', sep=''))     
   }
   if(filter=='WISE_W3'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_W3_WISE.tab', sep=''))     
   }
   if(filter=='WISE_W4'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_W4_WISE.tab', sep=''))     
  }
  
  if(filter==100 | filter=='Herschel_100'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_P100_Herschel.tab', sep=''))     
  }
   if(filter==160 | filter=='Herschel_160'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_P160_Herschel.tab', sep=''))     
   }
   if(filter==250 | filter=='Herschel_250'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_S250_Herschel.tab', sep=''))     
   }
   if(filter==350 | filter=='Herschel_350'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_S350_Herschel.tab', sep=''))     
   }
   if(filter==500 | filter=='Herschel_500'){
      out<-read.table(paste(FourXP_PATH,'/FourXP/data/calibrators/filters/filt_S450_Herschel.tab', sep=''))     
  }
  
 
  colnames(out)<-c('index','Wavelength', 'Trans')
  return=out
}
