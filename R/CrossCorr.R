#' CrossCorr
#'
#' @description Cross correlate one spectrum with a template
# "template" and "spec" must be the same length (even number). 
# Return a reduced size of array centred around zero shift. 
#' @param template 1D template spectrum 
#' @param spec 1D input spectrum 
#' @param corrSize size of output correlation vector
#' @param plan FFT plan
#' @author I. baldry, L. Davies <luke.j.davies@uwa.edu>
#' @examples 
#' None applicable as internal function.... 
#' @export
CrossCorr = function(template, spec, plan=plan, corrSize = 14000){

tempLength = length(template)

# take conjugate of fourier transform of rebinned template spectra
#fftTemplate <- FFT(template, plan = plan) / tempLength#TODO consider other (faster) transforms
fftTemplate <- fft(template)/ tempLength
fftTemplate <- Conj(fftTemplate)

#take fourier transform of spectra
#fftSpec <- FFT(spec, plan = plan) / tempLength
fftSpec <- fft(spec) / tempLength

# multiply by conj. of fourier transform of t_plate spectrum
fftSpec <- fftSpec * fftTemplate
#invfft <- IFFT(fftSpec, plan = plan, scale = FALSE)
invfft <- fft(fftSpec,inverse=TRUE)

#take real part of inverse FFT for cross correlation. 
crC <- Re(invfft) 
length <- length(crC)
halfLength <- length/2

crC <- c(crC[(halfLength+1):(length)], crC[1:(halfLength)])

#create output array: length = 2*corr_size+1
crossCorr <- crC[(halfLength-corrSize+1):(halfLength+corrSize+1)]
return = crossCorr
}
