#include <Rcpp.h>
using namespace Rcpp;

//This function returns the sign of a given real valued double.
// [[Rcpp::export]]
double signDblCPP (double x){
  double ret = 0;
  if(x > 0){ret = 1;}
  if(x < 0){ret = -1;}
  return(ret);
}

//Tested to be 6x faster(37 us vs 207 us). This operation is done from 200x per layer
//Original R function by Stas_G
// [[Rcpp::export]]
NumericVector findPeaksCPP( NumericVector vY, int m = 3) {
  int sze = vY.size();
  int i = 0;//generic iterator
  int q = 0;//second generic iterator
  
  int lb = 0;//left bound
  int rb = 0;//right bound
  
  bool isGreatest = true;//flag to state whether current index is greatest known value
  
  NumericVector ret(1);
  int pksFound = 0;
  
  for(i = 0; i < (sze-2); ++i){
    //Find all regions with negative laplacian between neighbors
    //following expression is identical to diff(sign(diff(xV, na.pad = FALSE)))
    if(signDblCPP( vY(i + 2)  - vY( i + 1 ) ) - signDblCPP( vY( i + 1 )  - vY( i ) ) < 0){
      //Now assess all regions with negative laplacian between neighbors...
      lb = i - m - 1;// define left bound of vector
      if(lb < 0){lb = 0;}//ensure our neighbor comparison is bounded by vector length
      rb = i + m + 1;// define right bound of vector
      if(rb >= (sze-2)){rb = (sze-3);}//ensure our neighbor comparison is bounded by vector length
      //Scan through loop and ensure that the neighbors are smaller in magnitude
      for(q = lb; q < rb; ++q){
        if(vY(q) > vY(i+1)){ isGreatest = false; }
      }
      
      //We have found a peak by our criterion
      if(isGreatest){
        if(pksFound > 0){//Check vector size.
          ret.insert( 0, double(i + 2) );
        }else{
          ret(0) = double(i + 2);
        }
        pksFound = pksFound + 1;
      }else{ // we did not find a peak, reset location is peak max flag.
        isGreatest = true;
      }//End if found peak
    }//End if laplace condition
  }//End loop
  return(ret);
}//End Fn