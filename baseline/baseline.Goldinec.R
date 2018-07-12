if (! require(foreach)) {
  install.packages('foreach')
  library(foreach)
}

T_rate <- function(x) {
  0.7679 + 11.2358*x - 39.7064*x^2 + 92.3583*x^3;
}

legend_c <- function(y, wl, p, s) {
  N <- length(wl)
  i <- order(wl)
  wl <- wl[i]
  y <- y[i]
  
  maxy <- max(y)
  dely <- ( maxy - min(y) ) / 2
  wl <- 2*(wl - wl[N]) / (wl[N] - wl[1]) + 1
  y <- (y - maxy) / dely + 1
  Tmat <- foreach(ord = 0:p, .combine = cbind) %do% {
    wl^ord
  }
  Tinv <- tcrossprod(MASS::ginv(crossprod(Tmat,Tmat)), Tmat)
  a <- Tinv %*% matrix(y,ncol = 1)
  z <- Tmat %*% a
  alpha <- 0.99*1/2          # Scale parameter alpha
  it    <- 0                 # Iteration number
  zp    <- matrix(1, N, 1)   # Previous estimation
  d     <- matrix(0, N, 1)
  while ( sum((z-zp)^2) / sum(zp^2) > 1e-9 ) {
    it  <- it+1   # Iteration number
    zp  <- z      # Previous estimation
    res <- y - z  # Residual
    for (num in 1:N){
      if (res[num] < s) {
        d[num] <- res[num]*(2*alpha-1)
      } else if (res[num] >= s) {
        d[num] <- -res[num] - alpha*(s^3) / (2*res[num]*res[num])
      }
    }
    a <- Tinv %*% (y + d) # Polynomial coefficients a
    z <- Tmat %*% a
  }
  j <- order(i)
  return( (z[j] - 1)*dely + maxy )
}

##' Goldindec semiautomated baseline correction method based on regression with cost fucntion.
##' 
##' @author Rustam Guliev
##' @param spc - spectra data in one of formats: data.frame, matrix, numeric vector or a hyperSpec object. If data.frame or matrix - one row is one spectrum. 
##' @param p - polynom order
##' @param peak_ratio - ratio of peaks
##' @param wl - optionally given wavelength values
##' @param eps - the parameter to terminate the iteration
##' @references Juntao Liu, Jianyang Sun, Xiuzhen Huang, Guojun Li, Binqiang Liu. Goldindec: A Novel Algorithm for Raman Spectrum Baseline Correction. DOI: 10.1366/14-07798

##' 
##' @export 
setGeneric ("baseline.Goldindec", function (spc, ...) standardGeneric ("baseline.Goldindec"))

##' @export
##' @rdname baseline.Goldindec
setMethod ("baseline.Goldindec", signature (spc = "numeric"), 
           function(spc, wl, p=4, peak_ratio=0.5, eps=0.0001) {
             stopifnot(length(spc) == length(wl), is.numeric(wl))
             udr <- T_rate(peak_ratio)
             
             a <- 0
             b <- 1
             s_old <- 1
             s <- a + 0.382*(b - a)
             z <- legend_c(spc, wl, p, s)
             up_down_rate <- sum(spc >= z) / sum(spc < z)             
             
             while ( (abs(up_down_rate - udr) > eps) && (abs(s_old - s) >= 0.00001) ) {
               s_old <- s
               if (up_down_rate - udr > eps) {
                 a <- s
               } else {
                 b <- s
               }
               s <- a + 0.382*(b - a)
               z <- legend_c(spc, wl, p, s)
               up_down_rate <- sum(spc >= z) / sum(spc < z)
             }
             
             return(z)
           })

##' @export
##' @rdname baseline.Goldindec
setMethod ("baseline.Goldindec", signature (spc = "matrix"),
           function(spc, wl = NULL,  ...) {
             if (is.null(wl)){
               wl <- hyperSpec::guess.wavelength(colnames(spc))
             }
             if (is.null(wl)) {
               stop("No wavelength data")
             }
             res <- apply(spc, 1, baseline.Goldindec, wl = wl, ...)
             if (nrow(res) != nrow(spc)){
               res <- t(res)
             }
             colnames(res) <- colnames(spc)
             rownames(res) <- rownames(spc)
             return(res)
           })

##' @export
##' @rdname baseline.Goldindec
setMethod ("baseline.Goldindec", signature (spc = "hyperSpec"),
           function(spc, ...) {
             validObject (spc)
             
             spc@data$spc <- baseline.Goldindec(unclass(spc@data$spc), wl = wl(spc), ...)
             return(spc)
           })
