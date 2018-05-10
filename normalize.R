#' Normalize spectra using different type of normalizations
#'
#' @param x - spectra data in one of formats: data.frame, matrix, numeric vector or a hyperSpec object. If data.frame or matrix - one row is one spectrum. 
#' @param type - type of normalization. 
#' `area` - makes area of each spectrum to be equal 1. 
#' `vector` - vector normalization, makes vector norm (i.e. sqrt(sum(\code{x}^2))) of each spectrum to be equal 1 
#' `minmax` - works the same as hyperSpec::normalize01: the values of a spectrum are mapped to [0, 1] by subtracting the minimum and subsequently dividing by the maximum. If all elements of \code{x} are equal, 1 is returned.
#' `max` - devides all values of a spectrum to max, i.e. makes maximum value to be qual 1.
#' `peak` - devides all values of a spectrum to a peak value, i.e. makes value of a peak to be qual 1. The value of peak is calculated as max value within a range of the peak defined by peak.range argument.  
#' @param peak.range - vector of two values defining range of peak, i.e. like \code{c(1500, 1600)}.
#' @param peak.index - logical. By default values in peak.range are related to colnames for data.frame and matrix, names for a vector, and wavelength for a hyperSpec object. However, indexes can be use by setting this argument to TRUE. 
#' @param na.rm - logical. Should missing values (including NaN) be removed?
#' @param tolerance tolerance level for determining what is 0 and 1 for 'minmax' normalization
#' @return Returns the normalized values in the format used in \code{x}, i.e. data.frame, matrix, hyperSpec or numeric vector. 
#' @author Rustam Guliev <glvrst@gmail.com>
#' @examples 
#' # Load test data from hyperSpec
#' library(hyperSpec)
#' flu.normalized <- normalize(flu, 'area')
#' flu.normalized <- normalize(flu, 'area')
#' laser.normalized <- normalize(laser, 'peak', peak.range=c(404,405.15))

normalize <- function (x, type = c('area', 'vector', 'minmax', 'max', 'peak'), 
                       peak.range = NULL, peak.index = FALSE, 
                       na.rm = TRUE, tolerance = hy.getOption ("tolerance")) {
  # Check class of X and make it a matrix
  if (is.data.frame(x)) {
    spc <- as.matrix(x)
  } else if (is.vector(x) && is.numeric(x)) {
    spc <- matrix(x, nrow = 1)
    colnames(spc) <- names(x)
  } else if (is (x,"hyperSpec")) {
    spc <- x$spc
    colnames(spc) <- wl(x)
  } else if (is.matrix(x)) {
    spc <- x
  } else {
    stop('Icorrect format of data. X must be one of following: "data.frame", "matrix", "numeric vector" or a "hyperSpec" object')
  }
  
  # Check type is correct. And peak.range is set if type is 'peak'
  type <- match.arg(type)
  
  switch(
    type,
    area = {
      m <- apply(spc, 1, mean, na.rm = na.rm)
      new.spc <- sweep (spc, 1, m, "/")
    },
    vector = {
      norm.to <- apply(spc, 1, function (x) {sqrt(sum(x*x, na.rm = na.rm))})
      new.spc <- sweep (spc, 1, norm.to, "/")
    },
    minmax = {
      m <- apply(spc, 1, min, na.rm = na.rm)
      new.spc <- sweep (spc, 1, m, `-`)
      m <- apply(new.spc, 1, max, na.rm = na.rm)
      new.spc <- sweep (new.spc, 1, m, `/`)
      new.spc [m < tolerance, ] <- 1
    },
    max = {
      m <- apply(spc, 1, function(x) {max(abs(x),na.rm = na.rm)})
      new.spc <- sweep (spc, 1, m, `/`)
    },
    peak = {
      stopifnot(is.vector(peak.range) && length(peak.range) == 2)
      if (peak.index) {
        pkrange <- min(peak.range):max(peak.range)
      } else {
        clnames <- as.numeric(colnames(spc))
        pkrange <- clnames >= min(peak.range) & clnames <= max(peak.range) 
      }
      norm.to <- apply(spc[, pkrange], 1, max, na.rm = na.rm)
      new.spc <- sweep (spc, 1, norm.to, "/")
    },
    {
      stop('Unknown type of normalization!')
    }
  )
  
  # Make output to have the same type as input x.
  if (is.data.frame(x)) {
    return(as.data.frame(new.spc))
  } else if (is.vector(x) && is.numeric(x)) {
    return(new.spc[1,])
  } else if (is (x,"hyperSpec")) {
    x$spc <- new.spc
    return(x)
  }
  return(new.spc) 
}