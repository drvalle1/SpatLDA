// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function makes multinomial draws
//Got this from https://stackoverflow.com/questions/24618370/using-rmultinom-with-rcpp
// [[Rcpp::export]]
IntegerVector rmultinom_1(NumericVector probs, int size) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(size, probs.begin(), k, ans.begin());
  return(ans);
}

//' This function samples zs when ys=1
// [[Rcpp::export]]
List SampleClustID(NumericMatrix theta, NumericMatrix phi, 
                   int nplot, int nspp, int ngrid, int nclust,
                   IntegerVector ArrayGSK,
                   IntegerVector ArraySoma,
                   IntegerMatrix ArrayTeste) {
  
  //convert array into arma::cube
  IntegerVector vecArray=clone(ArrayGSK);
  arma::cube ArrayGSK1(vecArray.begin(),ngrid, nspp, nclust);
  Integer Vector vecArray1=clone(ArraySoma);
  arma::cube ArraySoma1(vecArray1.begin(),ngrid, nspp, nplot);
  
  int n=0;
  NumericVector prob(ncommun);
  IntegerVector z(ncommun);
  
  for(int i=0; i<nplot; i++){
    for (int j=0; j<nspp; j++){
      n=ArrayTeste[j,i];
      if (n>0){
        for (int k=0; k<nclust; k++){
          prob(k)=theta(i,k)*phi(k,j);
        }
        prob=prob/sum(prob);
        for (int g=0; g<ngrid; g++){
          n=ArraySoma1[g,j,i];
          if (n>0){
            z=rmultinom_1(probs=prob, size=n);
            for (int k=0; k<nclust; k++){
              ArrayGSK(g,j,k)=z[k];
            }
          }
        }
      }
    }
  }
  return List::create(Named("ArrayGSK") = ArrayGSK);
}
