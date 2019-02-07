# TODO: Use R version if Rcpp is not available
Rcpp::sourceCpp('baseline/mor.cpp')

##' An automated baseline correction method based on iterative morphological operations.
##' 
##' @param spc - spectra data in one of formats: data.frame, matrix, numeric vector or a hyperSpec object. If data.frame or matrix - one row is one spectrum. 
##' @param w - half-width of moving window
##' @references Liankui Dai and YunliangChen. EXPRESS: An automated baseline correction method based on iterative morphological operations. DOI: 10.1177/0003702817752371
##' @author Rustam Guliev <glvrst@gmail.com>
##' 
##' @export 
setGeneric ("baseline.mor", function (spc, w) standardGeneric ("baseline.mor"))

##' @export
##' @rdname baseline.mor
setMethod ("baseline.mor", signature (spc = "numeric", w="numeric"), 
           function(spc, w) {
             stopifnot(length(w) == 1)
             res <- baseline_mor_cpp(spc, w)
             names(res) <- names(spc)
             return(res)
           })

##' @export
##' @rdname baseline.mor
setMethod ("baseline.mor", signature (spc = "matrix", w="numeric"), 
           function(spc, w) {
             stopifnot(length(w) == 1)
             res <- apply(spc, 1, baseline_mor_cpp, w = w)
             if (nrow(res) != nrow(spc)){
               res <- t(res)
             }
             colnames(res) <- colnames(spc)
             rownames(res) <- rownames(spc)
             return(res)
           })

##' @export
##' @rdname baseline.mor
setMethod ("baseline.mor", signature (spc = "hyperSpec", w="numeric"), 
           function(spc, w) {
             validObject (spc)
             stopifnot(length(w) == 1)
             spc@data$spc <- baseline.mor(unclass(spc@data$spc), w)
             return(spc)
           })

#---------------- I-MOR -----------------
##' @export 
setGeneric ("baseline.imor", function (spc) standardGeneric ("baseline.imor"))

##' @export
##' @rdname baseline.imor
setMethod ("baseline.imor", signature (spc = "numeric"), 
           function(spc) {
             res <- baseline_imor_cpp(spc)
             names(res) <- names(spc)
             return(res)
           })

##' @export
##' @rdname baseline.imor
setMethod ("baseline.imor", signature (spc = "matrix"), 
           function(spc) {
             res <- apply(spc, 1, baseline_imor_cpp)
             if (nrow(res) != nrow(spc)){
               res <- t(res)
             }
             colnames(res) <- colnames(spc)
             rownames(res) <- rownames(spc)
             return(res)
           })

##' @export
##' @rdname baseline.imor
setMethod ("baseline.imor", signature (spc = "hyperSpec"), 
           function(spc) {
             validObject (spc)
             spc@data$spc <- baseline.imor(unclass(spc@data$spc))
             return(spc)
           })

