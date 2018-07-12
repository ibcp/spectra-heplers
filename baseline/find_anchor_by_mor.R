source('baseline/baseline.mor.R')

##' Find anchor points of a spectrum using MOR-baseline low-variation ranges
##' 
##' An idea is to run MOR-baseline \code{n} times with \code{w} parameter varying from 1 to \code{n}. It's noticed that ranges around typical anchor points 
##' have a low variation of MOR-baseline. So, every range having a standard deviation of MOR-baseline lower than given \code{threshold} is expected to be a range around an anchor point.
##' Then anchor points themselves defined as local minimums within located ranges. 
##' 
##' @author Rustam Guliev
##' @param spc - numeric vector of a spectrum values
##' @param wl - numeric vector of wavelengths of the spectrum
##' @param n - number of iterations to run MOR-baseline
##' @param threshold - optionally given wavelength values
##' @return A matrix with two columns: (x, y) - respectively x and y coordinates of anchor points 


find_anchor_by_mor <- function(spc, wl, n = 50, threshold = 10) {
  # spc and wl are numeric vectors with the same length
  stopifnot(is.vector(spc), is.numeric(spc), is.vector(wl), is.numeric(wl), length(spc) == length(wl))
  
  # Get n MOR-baselines with w varying from 1 to n
  bl <- sapply(1:n, function(w) {baseline.mor(spc, w)})
  if (nrow(bl) != n) {
    bl <- t(bl)
  }
  
  # find MOR-baseline low-variance ranges
  morsd <- apply(bl, 2, sd)
  indx <- which(morsd<threshold)
  x <- wl[indx]
  y <- spc[indx]
  
  # define anchor points as local minimum 
  diff_left  <- diff(y)[1:(length(indx)-2)]
  diff_right <- diff(y)[2:(length(indx)-1)]
  indx_anchor <- which(diff_left<0 & diff_right>0)+1 
  
  return(cbind(x = x[indx_anchor], y = y[indx_anchor]))
}
