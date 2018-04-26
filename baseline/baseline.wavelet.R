##' Baseline using wavelets.
##'
##' @references Z.M. Zhang, S. Chen, Y.Z. Liang, et al., An intelligent background-correction algorithm for highly fluorescent samples in Raman spectroscopy. Journal of Raman Spectroscopy 41 (6), 659 (2010).

library(devtools); 
install_github("zmzhang/baselineWavelet")

baseline.wavelet <- function(spc, scales=1:70, wavelet='mexh', gapTh=3, skip=2, SNR.Th=1, ridgeLength=5,lambda=1000, differences=1) {
    apply(spc, 1, function(x) {
        wCoefs        <- cwt(x, scales=scales, wavelet=wavelet)
        localMax      <- getLocalMaximumCWT(wCoefs)
        ridgeList     <- getRidge(localMax, gapTh=gapTh, skip=skip)
        majorPeakInfo <- identifyMajorPeaks(x, ridgeList, wCoefs, SNR.Th=SNR.Th, ridgeLength=ridgeLength)
        peakWidth     <- widthEstimationCWT(x, majorPeakInfo)                                                
        return(baselineCorrectionCWT(x,peakWidth, lambda=lambda, differences=differences))
    })
}

