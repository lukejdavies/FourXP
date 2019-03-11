# FourXP (R package)

## Synopsis

FourXP (4XP) is the working R package for the 4MOST Extragalactic analysis pipeline. This includes tools for simulating spectra, observing with the 4MOST exposure time calculator, redshifting and pretty plotting tools.

## Installation

### Getting R

First things first, you will probably want to install a recent version of R that lets you build packages from source. The advantage of choosing this route is you can then update bleeding edge versions directly from GitHub. If you rely on the pre-build binaries on CRAN you might be waiting much longer.

#### Mac OS X

Pretty simple, go here and install the most relevant **pkg** file: <https://cloud.r-project.org/bin/macosx/>

#### Windows

Pretty simple, go here and install the most relevant **exe** file: <https://cloud.r-project.org/bin/windows/base/>

#### Linux

Varies a bit by platform. All the info on binaries is here: <https://cloud.r-project.org/bin/linux/>

Debian: `sudo apt-get install r-base r-base-dev`

Fedora: `sudo yum install R`

Suse:   More of a pain, see here <https://cloud.r-project.org/bin/linux/suse/README.html>

Ubuntu: `sudo apt-get install r-base-dev`

If you have a poorly supported version of Linux (e.g. CentOS) you will need to install R from source with the development flags (this bit is important). You can read more here: <https://cloud.r-project.org/sources.html>

Now you have the development version of R installed (hopefully) I would also suggest you get yourself R-Studio. It is a very popular and well maintained R IDE that gives you a lot of helpful shortcuts to scripting and analysing with R. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/> where you almost certainly want the free Desktop version.

### Getting FourXP


#### Source Install

The easiest way to install 4XP is via a source install in R. For this you will need the R devtools package: 

```R
install.packages('devtools')
library(devtools)
```

You can then install directly from github:

```R
install_github("lukejdavies/FourXP")
library(FourXP)
```

Following this, there are a number of model and template libraries you need to have installed. There is an inbuilt function in FourXP which will download and install these for you. However, you require ~1Gb of space to do this:

```R
install4XP(downloadModels = TRUE, ModelDir = ".")
```

Here 'ModelDir' is where you would like the models to be saved. Certain functions in FourXP require you to point to this directory (makeSpec() and FourXP_Sim()).

If you re-install FourXP, you will not need to re-download the libraries, but will need to run install4XP() as:

```R
install4XP(downloadModels = FALSE)
```

this unpacks various compressed files in the install.

#### R Package Dependencies

The above (binary or source install) should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **FourXP**:

```R
install.packages(c('FITSio', 'astro', 'magicaxis', 'Cairo', 'fftw'))
install.packages('devtools')
library(devtools)
install_github("lukejdavies/FourXP")
install4XP(downloadModels = FALSE)
```

#### External Dependencies

Some of the functionality of FourXP requires an installed version of the 4MOST Exposure Time Calculator (ETC). Functions such as observeSpec4FS() and FourXP_Sim() require that the ETC is excecutable from the command line as '4FS_ETC', you will also need to point towards the ETC system model directory (systemModDir) in these functions. 

Currently the ETC is only availabe to 4MOST team members. For further details contact L. Davies <luke.j.davies@uwa.edu>

If you are within 4MOST, you can get the ETC code here: http://wiki.4most.eu/4most-facility-simulator#toc8. Download the code and follow the install instructions in the README. Please then ensure that the ETC is executable from your command line as '4FS_ETC'. 

#### Examples

Assuming this has all installed successfully, you should now be able to load **FourXP** within R with the usual:

```R
library(FourXP)
```

You can get help and documentation for all FourXP functions using ?*function names*, i.e.:
```R
?makeSpec
```

## Generating a Model Spectrum

Here let's generate a spectrum at redshift, z, 0.56 with 19.8mag in the VST r-band. The template shape will come from a galaxy with g-i colour=0.55, stellar mass=10^10.2 and star-formation rate=7Msun/yr: 

```R
spec<-makeSpec(id='TestSpectrum', z=0.56, mag=19.8, band='VST_r', col=0.55, mass=10.2,sfr=7,agn='F')
```
A simple plot of the generated spectrum can be shown like this:

```R
plotSpec.basic(spec)
```

We can also make a more complicated plot with all of the information we input to the model:

```R
plotSpec(spec)
```

The z, magnitude, id, etc are stored in list 'spec' that we generated using makeSpec(). To look at the structure we have made try:

```R
names(spec)
```
Some of these have not yet been populated, such as expMin and prob.

Now try the same with a broad-line AGN model:

```R
specAGN<-makeSpec(id='TestAGNSpectrum', z=0.76, mag=17.6, band='VST_r', agn='B')
plotSpec(specAGN)
```

Using our first model we can run the redshifting code to see if we can automatically get the correct redshift (we should be able to as this spectrum has no noise!):

```R
FourXP_ZOut<-FourXP_Z(spec, verbose=F, doHelio=F)
FourXP_ZOut$results
         Z    Z1_PROB   CC_SIGMA   TEMPLATE         Z2  CC_SIGMA2  TEMPLATE2  CC_SIGMA3  CC_SIGMA4 
 0.5599835  1.0000000 15.3964229 43.0000000  0.1962446  3.6376233 47.0000000  3.0522262  2.8932451 
```

You will have had an error about and a lack of error, but don't worry about that for now! You can also see that we obtained a redshift of z=0.5599835 with a probability of 1, which is great given our input spectrum had z=0.56! Within 4XP, there is also a test version (this is very much in flux) of an emission-line pattern matching redshifting code FourXP_Zpat(). here we assume that the photo-z is known to 0.2 precision:  

```R
zFit<-FourXP_Zpat(spec, z_prior=c(0.4,0.6), plotCorr=F, filter=F)
Measured Redshift =  0.5596646 
Redshift Precision (in vs measured) =  64.46018 km/s 
```

So this gets roughly the correct redshift as well! 

For the next stages of this example, you will need to have the 4MOST ETC code installed....

Now we can mock observe this spectrum with 4MOST:

```R
observeSpec4FSOut<-observeSpec4FS(spec, expMin=120, nSub=4, keepFITS=F)
magplot(observeSpec4FSOut$blueRawWave, observeSpec4FSOut$blueRawFlux, type='l', col='blue', xlim=c(3500,10000), ylim=c(0, max(c(observeSpec4FSOut$blueRawFlux,observeSpec4FSOut$greenRawFlux,observeSpec4FSOut$redRawFlux),na.rm=T)), xlab='Wavelength, Ang', ylab='Counts')
lines(observeSpec4FSOut$greenRawWave, observeSpec4FSOut$greenRawFlux, type='l', col='darkgreen')
lines(observeSpec4FSOut$redRawWave, observeSpec4FSOut$redRawFlux, type='l', col='red')
```

Here we are running a 120min total exposure with 4 sub-exposures (of 30min each) and keeping all of the other observing defaults (see ?observeSpec4FS for available options). We can then plot the outputs for each spectral arm. We chose not to save the FITS image outputs of the ETC - `keepFITS=FALSE`.   

This code only provides the 3 spectral arms of 4MOST separately, we now want to combine them to one spectrum. We do this with stitch4MOST(), which can take in the output of observeSpec4FSOut(): 

```R
specObs<-stitch4MOST(observeSpec4FSOut=observeSpec4FSOut)
plotSpec.basic(specObs)
```

Or you can do this manually, saving the outputs of observeSpec4FS() as FITS files:

```R
observeSpec4FSOut<-observeSpec4FS(spec, expMin=120, nSub=4, keepFITS=T)
```
This makes a folder with *spec$id*_ETCout/ containing all of the 4FS_ETC outputs, which you can then point stitch4MOST() at:

```R
specObs<-stitch4MOST(id = 'TestSpectrum', blue_file = 'TestSpectrum_ETCout/specout_template_TestSpectrum_LRS_blue.fits', green_file = "TestSpectrum_ETCout/specout_template_TestSpectrum_LRS_green.fits", red_file = 'TestSpectrum_ETCout/specout_template_TestSpectrum_LRS_red.fits', blue_thru_file = '/Applications/4FS-ETC_app/4FS_ETC_system_model_v0.2/LRS/lrs_blue_material_4fs_efficiency_total.fits', green_thru_file = '/Applications/4FS-ETC_app/4FS_ETC_system_model_v0.2/LRS/lrs_green_material_4fs_efficiency_total.fits', red_thru_file = '/Applications/4FS-ETC_app/4FS_ETC_system_model_v0.2/LRS/lrs_red_material_4fs_efficiency_total.fits', ParaFile = "TestSpectrum_ETCout/ETC_input_params_tmp.txt", makePlot = F)
``

Finally, we can one again redshift our output spectrum to see if we would have obtained a redshift with this observation: 

```R
FourXP_ZOut<-FourXP_Z(specObs, verbose=F, doHelio=F)
FourXP_ZOut$results
          Z     Z1_PROB    CC_SIGMA    TEMPLATE          Z2   CC_SIGMA2   TEMPLATE2   CC_SIGMA3   CC_SIGMA4 
0.56019902  0.99967212 10.72456881 43.00000000  0.07671311  4.48136762 47.00000000  4.39744720  3.90459303 

FourXP_Zpat(specObs, z_prior=c(0.4,0.6), plotCorr=F, filter=T)
Measured Redshift =  0.5598869 
Redshift Precision (in vs measured) =  21.73643 km/s
```

And we got the right redshift back, so the observation would have been sucessful! Here we have set `filter=T` in `FourXP_Zpat()`. This is used to remove large noise spikes in the data prior to identifying emission lines (this wasn't needed in the earlier case, as the spectrum had no noise!).

All of this functionality is combined in the high level FourXP_Sim() function, which applied all stages of this process. You can input the desired spectrum you wish to produce and the observing contraints, and get the final observed spectrum back (including a lot of diagnostics!). Here is a full line example: 

```R
Out<-FourXP_Sim(id='Test', zIn=0.234, mag=19.8, band='VST_r', col=0.55, mass=10.2,sfr=7,agn='F', specDir='/Users/luke/work/IWG8/FourXPmodels/', expMin=60, nSub=3, SKYBRIGHT_TYPE='ZENITH', AIRMASS=1.4, IQ=1.1, SKYBRIGHT=21.77, TILT=6.0, MISALIGNMENT=0.1, systemModDir='/Applications/4FS-ETC_app/4FS_ETC_system_model_v0.2/', plot=T, verbose=TRUE)

Running FourXP_Sim, please wait..... 
     - Making simulated spectrum... 
     - Observing spectrum using 4FS ETC... 
     - Stitching Spectral Arms... 
     - Running 4XP_Z... 

 KEY INPUT PROPERTIES: 
 
   ID =  Test 
   Input Redshift =  0.234 
   ABMag =  19.8 
   MagBand =  VST_r 
   Template g-i col =  0.55 
   Template Log[M*] =  10.2 
   Template SFR =  7 
   Exposure Time (min) =  60 
   Number of sub-exposures =  3 
   Airmass =  1.4 
   Sky Brightness =  21.77 

 KEY OUTPUT PROPERTIES: 
 
Measured Redshift =  0.2340138 
Measured Probability =  0.9999997 
Redshift Precision (in vs measured) =  3.341352 km/s 
Signal to Noise Blue (median 4200-5000) =  3.078351 
Signal to Noise Green (median 5800-6600) =  5.370941 
Signal to Noise Red (median 7800-8600) =  6.844171 

 
FourXP_Sim finished! 
```

## Contributors

L. J. M. Davies, I. Baldry, L. Drygala, A. Robotham

## License

LGPL-3
