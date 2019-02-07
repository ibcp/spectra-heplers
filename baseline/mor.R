##' An automated baseline correction method based on iterative morphological operations.
##' 
##' @author Rustam Guliev
##' @param spc - matrix of spectra data: one row - one spectrum.
##' @param w - half-width of moving window
##' @references Liankui Dai and YunliangChen. EXPRESS: An automated baseline correction method based on iterative morphological operations. DOI: 10.1177/0003702817752371

# ============== SIMPLE Mor ==================
mor.erosion <- function(spc, w) {
  if (length(w) == 1) {
    w <- rep(w,nrow(spc))
  }
  newspc <- spc
  p <- ncol(spc)
  for (i in 1:nrow(spc)) {
    for (j in 1:p) {
      newspc[i,j] <- min(spc[i,max(1,j-w[i]):min(p,j+w[i])],na.rm = T)
    }
  }
  return(newspc)
}

mor.dilation <- function(spc, w) {
  if (length(w) == 1) {
    w <- rep(w,nrow(spc))
  }
  newspc <- spc
  p <- ncol(spc)
  for (i in 1:nrow(spc)) {
    for (j in 1:p) {
      newspc[i,j] <- max(spc[i,max(1,j-w[i]):min(p,j+w[i])],na.rm = T)
    }
  }
  return(newspc)
}

mor.opening <- function(spc,w){
  return(mor.dilation(mor.erosion(spc,w),w))
}

mor.P <- function(spc,w) {
  return(0.5*(mor.erosion(mor.opening(spc,w),w) + mor.dilation(mor.opening(spc,w),w)))
}

baseline.mor <- function(spc,w) {
  return(pmin(mor.P(spc,w),spc))
}

# =============== ITERATIVE Mor =============================
# Adaptively determine the optimal structuring element size
mor.getwopt <- function(spc,wstart=20) {
  wopts <- rep(NA,nrow(spc))
  for (i in 1:nrow(spc)) {
    s <- spc[i,,drop=FALSE]
    w <- wstart
    O <- mor.opening(s,w)
    diff <- Inf
    while(diff > 0.1) {
      prevO <- O
      w <- w+1
      O <- mor.opening(s,w)
      diff <- sum(abs(O-prevO))
    }
    wopts[i] <- w
  }
  return(wopts)
}

# Estimate the baseline based on iterative morphological operations
baseline.imor <- function(spc, tol=0.0001) {
  imor <- spc
  wopt <- mor.getwopt(spc)
  for (i in 1:nrow(spc)) {
    s <- spc[i,,drop=F]
    P <- mor.P(s, wopt[i])
    b <- pmin(P, s)
    RD <- sum((b-s)*(b-s))/sum(s*s)
    while(RD > tol) {
      newb <- pmin(mor.P(b, wopt[i]), s)
      RD <- sum((newb-b)*(newb-b))/sum(b*b)
      b <- newb
    }
    imor[i,] <- b
  }
  return(imor)
}
