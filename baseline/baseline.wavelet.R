if (! require(baselineWavelet)) {
  devtools::install_github("zmzhang/baselineWavelet")
  library(baselineWavelet)
}

##' Baseline using wavelets.
##'
##' @title Baseline using wavelets.
##' @param spc - spectra data in one of formats: data.frame, matrix, numeric vector or a hyperSpec object. If data.frame or matrix - one row is one spectrum. 
##' @return Returns baselines in the format used in \code{spc}, i.e. data.frame, matrix, hyperSpec or numeric vector. 
##' 
##' @references Z.M. Zhang, S. Chen, Y.Z. Liang, et al., An intelligent background-correction algorithm for highly fluorescent samples in Raman spectroscopy. Journal of Raman Spectroscopy 41 (6), 659 (2010).
##' @references GitHub: zmzhang/baselineWavelet
##' @author Rustam Guliev <glvrst@gmail.com>
##' 
##' @export 
setGeneric ("baseline.wavelet", function (spc, ...) standardGeneric ("baseline.wavelet"))

##' @export
##' @rdname baseline.wavelet
setMethod ("baseline.wavelet", signature (spc = "numeric"), 
           function(spc, scales=1:70, wavelet='mexh', gapTh=3, skip=2, SNR.Th=1, ridgeLength=5,lambda=1000, differences=1) {
             wCoefs        <- cwt(spc, scales=scales, wavelet=wavelet)
             localMax      <- getLocalMaximumCWT(wCoefs)
             ridgeList     <- getRidge(localMax, gapTh=gapTh, skip=skip)
             majorPeakInfo <- identifyMajorPeaks(spc, ridgeList, wCoefs, SNR.Th=SNR.Th, ridgeLength=ridgeLength)
             peakWidth     <- widthEstimationCWT(spc, majorPeakInfo)                                                
             return(baselineCorrectionCWT(spc, peakWidth, lambda=lambda, differences=differences))
           })

##' @export
##' @rdname baseline.wavelet
setMethod ("baseline.wavelet", signature (spc = "matrix"), 
           function(spc, ...) {
             res <- apply(spc, 1, baseline.wavelet, ... = ...)
             if (nrow(res) != nrow(spc)){
               res <- t(res)
             }
             colnames(res) <- colnames(spc)
             rownames(res) <- rownames(spc)
             return(res)
           })

##' @export
##' @rdname baseline.wavelet
setMethod (baseline.wavelet, signature (spc = "hyperSpec"), 
           function (spc, ...){
             validObject (spc)
             spc@data$spc <- baseline.wavelet(unclass(spc@data$spc), ...)
             return(spc)
           })  
