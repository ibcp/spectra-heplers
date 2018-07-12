##' A baseline correction method based splines on anchor points of a spectrum
##' 
##' @param x - whether a vector of x coordinates of anchor points, or matrix of two columns (x,y) - x and y coordinates of anchor points
##' @param y - a numeric vector of y coordinates. 
##' @param wl - a wavelength range of baseline
##' @return Returns a spline object with x and y coordinates of baseline.
##' @author Rustam Guliev <glvrst@gmail.com>
##' @example 
##' 
##' wl <- as.numeric(colnames(testspc))
##' anchors <- find_anchor_by_mor(testspc[1,], wl, n=50, threshold=25)
##' bl <- baseline.anchor(anchors[,1], anchors[,2], wl)
##' plot(x = wl, y=testspc[1,], type='l')
##' points(x = anchors[,1], y=anchors[,2], col='green',pch=15)
##' lines(x = bl$x, y=bl$y, col='red')

setGeneric ("baseline.anchor", function (x, y, wl) standardGeneric ("baseline.anchor"))

##' @export
##' @rdname baseline.anchor
setMethod ("baseline.anchor", signature (x = "numeric", y="numeric", wl="numeric"), 
           function(x, y, wl) {
             stopifnot( length(x) == length(y), is.vector(x), is.vector(y), is.vector(wl))
             stopifnot( any( wl>=min(x) & wl<=max(x) ) )
             if (min(wl)<min(x) | max(wl)>max(x)) {
               warning('There are some points in wl outside the range [min(x), max(x)]. wl range will be cutted.')
             }
             
             spline(x = x, y = y, xout = wl[wl>=min(x) & wl<=max(x)], method = "fmm")
           })