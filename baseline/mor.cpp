#include <Rcpp.h>
using namespace Rcpp;

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

// This function is simpy moving window min
NumericVector mor_erosion(const NumericVector x, const int w) {
  int n = x.size();
  NumericVector newx(n);
  unsigned long int w_start, w_end;
  
  
  for (int i = 0; i < n; i++) {
    if (NumericVector::is_na(x[i])) {
      newx[i] = NA_REAL;
    } else {
      w_start = ((i - w) < 0) ? 0 : i-w;
      w_end   = ((i + w) >= n) ? n-1 : i + w;
      newx[i] = vector_min(x, w_start, w_end);
    }
  }
  
  return newx;
}

// This function is simpy moving window max
NumericVector mor_dilation(const NumericVector x, const int w) {
  int n = x.size();
  NumericVector newx(n);
  unsigned long int w_start, w_end;
  
  
  for (int i = 0; i < n; i++) {
    if (NumericVector::is_na(x[i])) {
      newx[i] = NA_REAL;
    } else {
      w_start = ((i - w) < 0) ? 0 : i-w;
      w_end   = ((i + w) >= n) ? n-1 : i + w;
      newx[i] = vector_max(x, w_start, w_end);
    }
  }
  
  return newx;
}

NumericVector mor_opening(const NumericVector x, const int w) {
  return mor_dilation(mor_erosion(x,w),w);
}

NumericVector mor_P(const NumericVector x, const int w) {
  return 0.5*(mor_erosion(mor_opening(x,w),w) + mor_dilation(mor_opening(x,w),w));
}

// Min of mor_P and original x
// [[Rcpp::export]]
NumericVector baseline_mor_cpp(const NumericVector x, const int w) {
  int n = x.size();
  NumericVector bl = mor_P(x, w);
  for(int i=0;i<n;i++){
    if (x[i] < bl[i]) {
      bl[i] = x[i];
    }
  }
  return bl;
}

// =============== ITERATIVE Mor =============================
// Adaptively determine the optimal structuring element size
// [[Rcpp::export]]
unsigned int mor_getwopt(NumericVector x, const unsigned int wstart=10) {
  
  unsigned int w = wstart;
  int n = x.size();
  NumericVector old_opening(n);
  NumericVector new_opening = mor_opening(x, w);
  double diff = 10000;
  
  while (diff > 0) {
    w++;
    old_opening = new_opening;
    new_opening = mor_opening(x, w);
    diff = sum(abs(new_opening-old_opening));
  }
  
  return w-1;
}

// Estimate the baseline based on iterative morphological operations
// [[Rcpp::export]]
NumericVector baseline_imor_cpp(const NumericVector x, const double tol=0.0001){
  unsigned int w = mor_getwopt(x);
  NumericVector b(x.size());
  std::copy( x.begin(), x.end(), b.begin() ); // b = x;
  NumericVector b_new = pmin(mor_P(b, w),x);
  double rd = sum((b_new-b)*(b_new-b)) / sum(b*b);
  
  while(rd > tol) {
    std::copy( b_new.begin(), b_new.end(), b.begin() ); // b = b_new;
    b_new = pmin(mor_P(b, w), x);
    rd = sum((b_new-b)*(b_new-b)) / sum(b*b);
  }
  return b_new;
}

