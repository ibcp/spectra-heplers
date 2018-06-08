#include <Rcpp.h>
// #include <algorithm>
using namespace Rcpp;

// TODO:
// - Handle NA's. Use std::isnan
// - Moving window min/max.
// - I-MOR
// - parallel
// - matrix

// [[Rcpp::export]]
double vector_min(const NumericVector x, const int start, const int end) {
  int i = 0;
  double smallest = x[start];
  
  while ((start+i) <= end){
    if (x[start+i] < smallest) {
      smallest = x[start+i];
    }
    i++;
  }
  
  return smallest;
}

// [[Rcpp::export]]
double vector_max(const NumericVector x, const int start, const int end) {
  int i = 0;
  double biggest = x[start];
  
  while ((start+i) <= end){
    if (x[start+i] > biggest) {
      biggest = x[start+i];
    }
    i++;
  }
  
  return biggest;
}

// [[Rcpp::export]]
NumericVector mor_erosion(const NumericVector x,
                                const int w) {
  int n = x.size();
  NumericVector newx(n);
  unsigned long int w_start, w_end;
  
  
  for (int i = 0; i < n; i++) {
    w_start = ((i - w) < 0) ? 0 : i-w;
    w_end   = ((i + w) >= n) ? n-1 : i + w;
    newx[i] = vector_min(x, w_start, w_end);
  }
  
  return newx;
}

// [[Rcpp::export]]
NumericVector mor_dilation(const NumericVector x, const int w) {
  int n = x.size();
  NumericVector newx(n);
  unsigned long int w_start, w_end;
  
  
  for (int i = 0; i < n; i++) {
    w_start = ((i - w) < 0) ? 0 : i-w;
    w_end   = ((i + w) > n) ? x.size() : i + w;
    newx[i] = vector_max(x, w_start, w_end);
  }
  
  return newx;
}

// [[Rcpp::export]]
NumericVector mor_opening(const NumericVector x, const int w) {
  return mor_dilation(mor_erosion(x,w),w);
}

// [[Rcpp::export]]
NumericVector mor_P(const NumericVector x, const int w) {
  return 0.5*(mor_erosion(mor_opening(x,w),w) + mor_dilation(mor_opening(x,w),w));
}

// [[Rcpp::export]]
NumericVector baseline_mor(const NumericVector x, const int w) {
  int n = x.size();
  NumericVector bl = mor_P(x, w);
  for(int i=0;i<n;i++){
    if (x[i] < bl[i]) {
      bl[i] = x[i];
    }
  }
  return bl;
}

/*
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
 */

