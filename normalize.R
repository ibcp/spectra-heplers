.area <- function(x, na.rm=TRUE) {
  wl <- hyperSpec::guess.wavelength(names(x))
  if (is.null(wl)) {
    warning('No wavelengths are provided. 1:lenght(x) is used instead.')
    wl <- 1:length(x)
  }
  if (na.rm) {
    clr <- !is.na(x)
    AUC <- sum(0.5 * diff(wl[clr]) * (head(x[clr],-1) + tail(x[clr],-1)))
  } else {
    if (anyNA(x)) {
      warning('NA has occured. All values are set to NA.')
      AUC <- NA
    } else {
      AUC <- sum(0.5 * diff(wl) * (head(x,-1) + tail(x,-1)))
    }
  }
  return(AUC)
}

#' Normalize spectra using different type of normalizations
#'
#' @title Normalize spectra
#' @param x - spectra data in one of formats: data.frame, matrix, numeric vector or a hyperSpec object. If data.frame or matrix - one row is one spectrum. 
#' @param type - type of normalization. 
#' `area` - makes area of each spectrum to be equal 1.
#' `mean` - makes mean of each spectrum to be equal 1. This may be used as an alternative for the area normalization since there is no big difference between these two types of normalization, but mean is faser to calculate.
#' `vector` - vector normalization, makes vector norm (i.e. sqrt(sum(\code{x}^2))) of each spectrum to be equal 1 
#' `zeroone` - works the same as hyperSpec::normalize01: the values of a spectrum are mapped to [0, 1] by subtracting the minimum and subsequently dividing by the maximum. If all elements of \code{x} are equal, 1 is returned.
#' `max` - devides all values of a spectrum to max, i.e. makes maximum value to be qual 1.
#' `peak` - devides all values of a spectrum to a peak value, i.e. makes value of a peak to be qual 1. The value of peak is calculated as max value within a range of the peak defined by peak.range argument.  
#' @param peak.range - vector of two values defining range of peak, i.e. like \code{c(1500, 1600)}.
#' @param peak.index - logical. By default values in peak.range are related to colnames for data.frame and matrix, names for a vector, and wavelength for a hyperSpec object. However, indexes can be use by setting this argument to TRUE. 
#' @param na.rm - logical. Should missing values (including NaN) be ingnored while calculating normalization coefficient, i.e. sum, mean, area, etc.?
#' @param tolerance tolerance level for determining what is 0 and 1 for 'zeroone' normalization
#'
#' @return Returns the normalized values in the format used in \code{x}, i.e. data.frame, matrix, hyperSpec or numeric vector. 
#' @author Rustam Guliev <glvrst@gmail.com>
#' 
#' @examples 
#' # Load test data from hyperSpec
#' library(hyperSpec)
#' flu.normalized <- normalize(flu, 'area')
#' flu.normalized <- normalize(flu, 'vector')
#' laser.normalized <- normalize(laser, 'peak', peak.range=c(404,405.15))
#' 
#' @export 
setGeneric ("normalize", function (x, ...) standardGeneric ("normalize"))

##' @export
##' @rdname normalize
setMethod ("normalize", signature (x = "matrix"), 
           function (x, 
                     type = c('area', 'mean', 'vector', 'zeroone', 'max', 'peak'), 
                     peak.range = NULL, peak.index = FALSE, 
                     na.rm = TRUE, tolerance = hy.getOption ("tolerance")) {
             # Check type is correct. And peak.range is set if type is 'peak'
             type <- match.arg(type)
             
             switch(
               type,
               area = {
                 m <- apply(x, 1, .area, na.rm = na.rm)
                 new.x <- sweep (x, 1, m, "/")
               },
               mean = {
                 m <- rowMeans(x, na.rm = na.rm)
                 new.x <- sweep (x, 1, m, "/")
               },
               vector = {
                 norm.to <- apply(x, 1, function (x) {sqrt(sum(x*x, na.rm = na.rm))})
                 new.x <- sweep (x, 1, norm.to, "/")
               },
               zeroone = {
                 m <- apply(x, 1, min, na.rm = na.rm)
                 new.x <- sweep (x, 1, m, `-`)
                 m <- apply(new.x, 1, max, na.rm = na.rm)
                 new.x <- sweep (new.x, 1, m, `/`)
                 new.x [m < tolerance, ] <- 1
               },
               max = {
                 m <- apply(x, 1, function(x) {max(abs(x),na.rm = na.rm)})
                 new.x <- sweep (x, 1, m, `/`)
               },
               peak = {
                 stopifnot(is.vector(peak.range) && length(peak.range) == 2)
                 if (peak.index) {
                   pkrange <- min(peak.range):max(peak.range)
                 } else {
                   clnames <- as.numeric(colnames(x))
                   pkrange <- clnames >= min(peak.range) & clnames <= max(peak.range) 
                 }
                 norm.to <- apply(x[, pkrange], 1, max, na.rm = na.rm)
                 new.x <- sweep (x, 1, norm.to, "/")
               },
               {
                 stop('Unknown type of normalization!')
               }
             )
             return(new.x)
           })

##' @export
##' @rdname normalize
setMethod (normalize, signature (x = "numeric"), function (x, ...){
  # Convert vector to a matrix of one row. 
  spc <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  res <- normalize(spc, ...)
  # Convert result matrix to a vector
  return(res[1,])
})

##' @export
##' @rdname normalize
setMethod (normalize, signature (x = "hyperSpec"), function (x, ...){
  validObject (x)
  
  colnames(x@data$spc) <- wl(x)
  x@data$spc <- normalize (unclass (x@data$spc), ...)
  
  ## logbook
  return(x)
})