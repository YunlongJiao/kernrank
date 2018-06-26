#include <Rcpp.h>
using namespace std;

// This is the core implementation of weighted Kendall kernel (kendall_weight_inner.R).
// It is essentially a modified quicksort algorithm that counts weighted non-inversion number of a permutation.

int kendall_weight_quickC(double *noninv, 
                          vector<int>& idxC, 
                          vector<int>& zC, 
                          int nC, 
                          int keyC, 
                          int kC, 
                          vector<double>& uC){
  
  if (idxC.size() > 1) {
    vector<int> idxChigh, idxClow;
    int ipiv = std::rand() % idxC.size();
    double zpiv = zC[idxC[ipiv] - 1], cnum = 0, cmin = 0, ctop = 0, cweight = 0, cconvweight = 0, cmultweight = 0;

    for (vector<int>::iterator it = idxC.begin(); it != idxC.end(); ++it) {
      if (zC[*it - 1] < zpiv) {
        idxClow.push_back(*it);
        cnum += 1.0;
        cmin += min(*it, zC[*it - 1]) / (double)nC;
        ctop += (((*it >= kC) && (zC[*it - 1] >= kC)) ? 1.0 : 0.0);
        cweight += uC[*it - 1];
        cconvweight += uC[zC[*it - 1] - 1];
        cmultweight += uC[*it - 1] * uC[zC[*it - 1] - 1];
      } else {
        idxChigh.push_back(*it);
        if (keyC == 1) {
          // "ken"
          *noninv += cnum;
        } else if (keyC == 2) {
          // "aken"
          *noninv += cmin;
        } else if (keyC == 3) {
          // "top"
          *noninv += (((*it >= kC) && (zC[*it - 1] >= kC)) ? ctop : 0.0);
        } else if (keyC == 4) {
          // "add"
          *noninv += cmultweight + cweight * uC[zC[*it - 1] - 1] + cconvweight * uC[*it - 1] + cnum * uC[*it - 1] * uC[zC[*it - 1] - 1];
        } else if (keyC == 5) {
          // "mult"
          *noninv += cmultweight * uC[*it - 1] * uC[zC[*it - 1] - 1];
        }
      }
    }
    
    // for (vector<int>::const_iterator i = idxChigh.begin(); i != idxChigh.end(); ++i) {
    //   Rcpp::Rcout << " high " << *i << endl;
    // }
    kendall_weight_quickC(noninv, idxChigh, zC, nC, keyC, kC, uC);
    // Rcpp::Rcout << " noninv " << *noninv << endl;
    
    // for (vector<int>::const_iterator i = idxClow.begin(); i != idxClow.end(); ++i) {
    //   Rcpp::Rcout << " low " << *i << endl;
    // }
    kendall_weight_quickC(noninv, idxClow, zC, nC, keyC, kC, uC);
    // Rcpp::Rcout << " noninv " << *noninv << endl;
  }
  
  return 0;
}

// [[Rcpp::export]]
double kendall_weight_quickR(Rcpp::IntegerVector idxR, 
                             Rcpp::IntegerVector zR, 
                             int nR, 
                             int keyR, 
                             int kR, 
                             Rcpp::NumericVector uR){
  double res = 0;
  vector<int> idxC = Rcpp::as<vector<int> >(idxR);
  vector<int> zC = Rcpp::as<vector<int> >(zR);
  int nC = nR, keyC = keyR, kC = kR;
  vector<double> uC = Rcpp::as<vector<double> >(uR);
  
  kendall_weight_quickC(&res, idxC, zC, nC, keyC, kC, uC);
  
  return res;
}
