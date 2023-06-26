// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;     

// CRPS divergence for uniform distribution
// [[Rcpp::export]]
double crps_div(arma::colvec x){
  
  double out1 = 0;
  double out2 = 0;
  double n = x.size();
  for (int i = 0; i < n; i++) {
    out1 += pow(x[i], 2) - x[i] + 0.5;
    for (int j = i; j < n; j++) {
      out2 += abs(x[i] - x[j]);
    }
  }
  double out = out1/n - out2/pow(n, 2) - 1.0/6.0;
  return(out);
}

