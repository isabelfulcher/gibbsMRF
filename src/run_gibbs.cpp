#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <random>
using namespace Rcpp;

float expit (float x){

  float output = exp(x)/(1+exp(x));
  return output;

}

float neighborsum (NumericVector v1, NumericVector v2, float weight, int N){

  int sum = 0;

  for (int i=0; i<N; i++) {
    sum += (v1[i]*v2[i]);
  }
  float scaled = sum/weight;

  return scaled;

}


float bernoulli (float p){
  float r = ((float) rand() / (RAND_MAX));
  unsigned int br = 0;
  if (r >= p)
    br = 1;
  return br;
}



//' @export
// [[Rcpp::export]]
IntegerMatrix run_gibbSimple(NumericMatrix adj, NumericVector weights, float alpha0, float alpha1, int R, int N, NumericVector start) {

  IntegerMatrix mat(N, R);
  NumericVector vec = start;

  for (int r=0; r<R; r++) {
    for (int i=0; i<N; i++) {
    float sum = neighborsum(vec,adj(i,_),weights(i),N);
    float prob = expit(alpha0 + alpha1*sum);
    vec[i] = bernoulli(prob);
    }
    mat(_,r) = vec ;
  }
  return mat;
}

