#' Create PCA-LDA model able to predict using original variables instead of PCs.
#'
#' @title PCA-LDA combination
#' @param x - original data to be transformed by PCA 
#' @param y - target (grouping) variable for LDA
#' @param PCs - integer defining number of principal components (PCs) OR integer vector of PCs  to use in LDA 
#' @param center - center parameter to be sent to prcomp
#' @param scale - scale. parameter to be sent to prcomp
#' @param prior - vector of prior probabilities to use in lda
#'
#' @return Returns a list/pclda object of the following components:
#'  PCs - principal components used for LDA
#'  center  - center used in PCA
#'  scale - scale used in PCA
#'  prior - priors used in LDA
#'  scaling - rotation matrix/vector to project newdata to the space of LD(s)
#'  means - group means in the space of original variables
#'  lev - levels of grouping (\code{y}) variable
#' @author Rustam Guliev <glvrst@gmail.com>
#' 
#' @examples 
#' pclda_model <- pclda(iris[,-5],iris[,5],PCs = 2)
#' pr <- predict.pclda(pclda_model, iris[,-5])
#' table(iris$Species, pr$class)
#' plot(pr$x[,1], pr$x[,2], pch=16, col=iris$Species, xlab='LD 1', ylab = 'LD 2')

pclda <- function(x, y, PCs, center=TRUE, scale=FALSE, prior=NULL) {
  # Prepare parameters
  if (length(PCs) == 1) {
    PCs <- 1:PCs
  }
  if (! is.factor(y)){
    y <- factor(y)
  }
  if (is.null(prior)) {
    prior <- rep(1/nlevels(y), nlevels(y))
    names(prior) <- levels(y)
  }
  
  # PCA + LDA
  pca <- prcomp(x, center = center, scale. = scale)
  X <- pca$x[,PCs,drop=FALSE]
  l <- lda(x = X, grouping=y, prior = prior)
  
  # Prepare information for output. This is needed for predictions
  means <- as.matrix(aggregate(x, by=list(y), mean, na.rm=TRUE)[,-1])
  rownames(means) <- levels(y)
  colnames(means) <- colnames(x)
  
  # Return
  res <- list(
    PCs     = PCs,
    center  = pca$center,
    scale   = pca$scale,
    prior   = l$prior,
    scaling = pca$rotation[,PCs] %*% l$scaling[,,drop=F],
    means   = means,
    lev     = levels(y)
  )
  class(res) <- append(class(res),"pclda")
  return(res)
}

predict.pclda <- function(object, newdata, prior = object$prior) {
  
  ng <- length(object$prior)
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1) 
      stop("invalid 'prior'")
    if (length(prior) != ng) 
      stop("'prior' is of incorrect length")
  }
  
  # PROJECT newdata to LDs
  x  <- scale(newdata, center = object$center, scale = object$scale) %*% object$scaling
  dm <- scale(object$means, center = object$center, scale = object$scale) %*% object$scaling
  
  # Calculate posteriors
  dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x), length(prior), byrow = TRUE) - x %*% t(dm)
  dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
  posterior <- dist/drop(dist %*% rep(1, ng))
  
  # Define predicted class
  nm <- names(object$prior)
  cl <- factor(nm[max.col(posterior)], levels = object$lev)
  dimnames(posterior) <- list(rownames(x), nm)
  
  list(class = cl, posterior = posterior, x = x)
}