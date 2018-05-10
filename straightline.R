#' Generate straight line between two points
#' 
#' The values in given data are replaced by straight line. 
#' This is one of ways that cosmic particles cotribution in Raman/SERS spectra can be handled.
#'
#' @param data - spectra data in one of formats: data.frame, matrix, numeric vector or a hyperSpec object. If data.frame or matrix - one row is one spectrum.
#' @param x - vector of x axis values. By defauld colnames/names used.
#' @param from - defines x value of left side of the range to be strainghtlined. 
#' @param to - defines x value of right side of the range to be strainghtlined.
#' @param index - logical. Should values \code{from} and \code{to} to be interpreted as indexes? 
#' @return Matrix with values replaced with a strainght line.
#' @author Rustam Guliev <glvrst@gmail.com>
#' @examples 
#' # Load test data from hyperSpec
#' library(hyperSpec)
#' 
#' # Pretend we need to replace range between 900 and 1200 for the first ten spectra.
#' res <- straightline(chondro[1:10], from = 990, to = 1010)
#' plotspc(res[,,900~1100])

straightline <- function(data, x = NULL, from = NULL, to = NULL, index = FALSE) {
  # Define x axis values
  if (is.null(x)) {
    if (is.data.frame(data) || is.matrix(data)) {
      x <- as.numeric(colnames(data))
    } else if (is.vector(data) && is.numeric(data)) {
      x <- as.numeric(names(data))
    } else if (is (data,"hyperSpec")) {
      x <- wl(data)
    } else {
      stop('Icorrect format of data. DATA must be one of following: "data.frame", "matrix", "numeric vector" or a "hyperSpec" object')
    }
  }
  
  # Define the range to be replaced by straight line
  # This is simply defining indexes of "left" and "right" points of the line.
  # cols - defines columns to be changed. Usualy they are just columns between 'left' and 'right' indexed,
  # but we make it for flexible in case columns are not ordered.
  if (index) {
    if (is.null(from)) {
      from <- 1
    }
    if (is.null(to)) {
      to <- length(x)
    }
    indx.left <- from
    indx.right <- to
    cols <- from:to
  } else {
    if (is.null(from)) {
      from <- min(x)
    }
    if (is.null(to)) {
      to <- max(x)
    }
    indx.left <- which(x == from)
    if (length(indx.left) == 0) {
      stop(paste0('No such wavelength "', from, '". Use existing wavelength for from and to parameters'))
    }
    indx.right <- which(x == to)
    if (length(indx.right) == 0) {
      stop(paste0('No such wavelength "', to, '". Use existing wavelength for from and to parameters'))
    }    
    cols <- (x >= from) & (x <= to)
  }
   
  # Use formula y = k*x + b to generate a straight line. 
  # Having two points (left and right) (x1,y1) and (x2,y2), k = (y2-y1)/(x2-x1) and b = (x2*y1 - x1*y2)/(x2-x1)
  x1 <- x[indx.left]
  x2 <- x[indx.right]
  if (is (data,"hyperSpec")) {
    y1 <- data$spc[,indx.left]
    y2 <- data$spc[,indx.right]
  } else if (is.vector(data)) {
   y1 <- data[indx.left]
   y2 <- data[indx.right]
  } else {
   y1 <- data[,indx.left]
   y2 <- data[,indx.right]
  }

  k <- matrix((y2-y1) / (x2-x1), ncol=1)
  b <- (x2*y1 - x1*y2) / (x2-x1)
  
  if (is (data,"hyperSpec")) {
    data$spc[,cols] <- k %*% x[cols] + b
  } else if (is.vector(data)) {
    data[cols] <- k[,1] * x[cols] + b
  } else {
    data[,cols] <- k %*% x[cols] + b
  }
  return(data)
}

