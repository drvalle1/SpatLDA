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
IntegerVector rmultinom1(NumericVector probs, int size) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(size, probs.begin(), k, ans.begin());
  return(ans);
}

//' This function samples the cluster id of each tree
// [[Rcpp::export]]
List SampleClustID(NumericMatrix theta, NumericVector phi, 
                   int nplot, int ngrid, int nclust,
                   NumericVector ArrayGKP) {
  
  //convert array into arma::cube
  NumericVector vecArray=clone(ArrayGKP);
  arma::cube ArrayGKP1(vecArray.begin(),ngrid, nclust,nplot);
  
  //get SomaGP and SomaP
  IntegerMatrix SomaGP(ngrid,nplot);
  IntegerVector SomaP(nplot);
  for(int i=0; i<nplot; i++){
    for (int g=0; g<ngrid; g++){
      for (int k=0; k<nclust; k++){
        SomaGP(g,i)=SomaGP(g,i)+ArrayGKP1(g,k,i);
        SomaP[i]=SomaP[i]+ArrayGKP1(g,k,i);
      }
    }
  }
  
  //sample ClustID
  int n=0;
  NumericVector prob(nclust);
  IntegerVector z(nclust);

  for(int i=0; i<nplot; i++){
    n=SomaP[i];
    if (n>0){
      for (int k=0; k<nclust; k++){
        prob(k)=theta(i,k)*phi(k);
      }
      prob=prob/sum(prob);
      for(int g=0; g<ngrid; g++){
        n=SomaGP(g,i);
        if (n>0){
          z=rmultinom1(prob,n);
          for (int k=0; k<nclust; k++){
            ArrayGKP1(g,k,i)=z[k];
          }
        }
      }
    }
  }
  return List::create(Named("ArrayGKP1") = ArrayGKP1);
}

//' This function samples the plot id of each tree
// [[Rcpp::export]]
List SamplePlotID(NumericMatrix theta, NumericMatrix delta, 
                  int ngrid, int nclust, int nplot,
                  NumericVector ArrayGKP) {
  
  //convert array into arma::cube
  NumericVector vecArray=clone(ArrayGKP);
  arma::cube ArrayGKP1(vecArray.begin(),ngrid, nclust,nplot);
  
  //get SomaGK and SomaK
  IntegerMatrix SomaGK(ngrid,nclust);
  IntegerVector SomaK(nclust);
  for(int i=0; i<nplot; i++){
    for (int g=0; g<ngrid; g++){
      for (int k=0; k<nclust; k++){
        SomaGK(g,k)=SomaGK(g,k)+ArrayGKP1(g,k,i);
        SomaK[k]=SomaK[k]+ArrayGKP1(g,k,i);
      }
    }
  }
  
  //sample PlotID
  int n=0;
  NumericVector prob(nplot);
  IntegerVector z(nplot);
  
  for(int k=0; k<nclust; k++){
    n=SomaK[k];
    if (n>0){
      for (int g=0; g<ngrid; g++){
        n=SomaGK(g,k);
        if (n>0){
          for (int i=0; i<nplot; i++){
            prob[i]=theta(i,k)*delta(g,i);
          }
          prob=prob/sum(prob);
          z=rmultinom1(prob,n);
          for (int i=0; i<nplot; i++){
            ArrayGKP1(g,k,i)=z[i];  
          }
        }
      }
    }
  }
  return List::create(Named("ArrayGKP1") = ArrayGKP1,
                      Named("SomaGK")=SomaGK);
}

//' This function calculates SumPK
// [[Rcpp::export]]
IntegerMatrix GetSumPK(int nplot, int nclust, int nspp,
                    NumericVector SomaSKP) {
  //convert array into arma::cube
  NumericVector vecArray=clone(SomaSKP);
  arma::cube SomaSKP1(vecArray.begin(),nspp, nclust,nplot);
  
  IntegerMatrix SomaPK(nplot,nclust);
  for(int i=0; i<nplot; i++){
    for (int s=0; s<nspp; s++){
      for (int k=0; k<nclust; k++){
        SomaPK(i,k)=SomaPK(i,k)+SomaSKP1(s,k,i);
      }
    }
  }
  
  return SomaPK;
}

//' This function calculates SumKS
// [[Rcpp::export]]
IntegerMatrix GetSumKS(int nplot, int nclust, int nspp,
                       NumericVector SomaSKP) {
  //convert array into arma::cube
  NumericVector vecArray=clone(SomaSKP);
  arma::cube SomaSKP1(vecArray.begin(),nspp, nclust,nplot);
  
  IntegerMatrix SomaKS(nclust,nspp);
  for(int i=0; i<nplot; i++){
    for (int s=0; s<nspp; s++){
      for (int k=0; k<nclust; k++){
        SomaKS(k,s)=SomaKS(k,s)+SomaSKP1(s,k,i);
      }
    }
  }
  
  return SomaKS;
}
