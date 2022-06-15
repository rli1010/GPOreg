#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double sumC(NumericVector xb,NumericVector pihat,NumericVector wi) {
  int n = xb.size();   double total = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if(xb[i]>=xb[j] && i!=j) total=total+pihat[i]*wi[i]*wi[j];
    }}
  total=total/(n*(n-1));
  return total;  }


// [[Rcpp::export]]
double kernelC(NumericVector xb,NumericVector pihat,NumericVector wi,double Cn) {
  int n = xb.size();   double total = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if(i!=j){
        double addi=0.5 * erfc(-(xb[i]-xb[j]) * M_SQRT1_2/Cn);
        total=total+pihat[i]*wi[i]*wi[j]*addi;
      }}}
  total=total/(n*(n-1));
  return total;  }


// [[Rcpp::export]]
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

// [[Rcpp::export]]
double kernelC2(NumericVector xb,NumericVector pihat,NumericVector wi,double Cn) {
  int n = xb.size();   double total = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      if(i!=j){
         double addi=phi((xb[i]-xb[j])/Cn);
         total=total+pihat[i]*wi[i]*wi[j]*addi+pihat[j]*wi[i]*wi[j]*(1-addi);
      }}}
  total=total/(n*(n-1));
  return total;  }

// [[Rcpp::export]]
NumericVector myrank(NumericVector x,NumericVector wi) {
  int n=x.size();   NumericVector out(n);
  
  for (int i = 0; i < n; i++) {
    out[i]=0.5;
    for (int j = 0; j < n; j++) {
      if(x[j]<=x[i]) out[i]=out[i]+0.5*wi[j];
      if(x[j]<x[i])  out[i]=out[i]+0.5*wi[j];
    }}
  return out; }