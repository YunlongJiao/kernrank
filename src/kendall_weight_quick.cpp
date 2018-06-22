#include <Rcpp.h>
using namespace Rcpp;

// This is the core implementation of weighted Kendall kernel (kendall_weight_inner.R).
// It is essentially a modified quicksort algorithm that counts weighted non-inversion number of a permutation.

// [[Rcpp::export]]
double kendall_weight_quick(IntegerVector idx, 
                            const IntegerVector zC, 
                            const int nC, 
                            const int keyC, 
                            const Nullable<int> kC, 
                            const Nullable<NumericVector> uC) {
  static double noninv = 0;
  
  if (idx.size() > 0) {
    std::vector<int> idxhigh, idxlow;
    double zpiv = zC[idx[0] - 1], cnum = 0, cmin = 0, ctop = 0, cweight = 0, cconvweight = 0, cmultweight = 0;
    
    for (IntegerVector::iterator it = idx.begin(); it != idx.end(); ++it) {
      if (zC[*it - 1] < zpiv) {
        idxlow.push_back(*it);
        cnum += 1.0;
        cmin += std::min(*it, zC[*it - 1]) / (double)nC;
        if (kC.isNotNull()) {
          ctop += (((*it >= kC) && (zC[*it - 1] >= kC)) ? 1.0 : 0.0);
        }
        if (uC.isNotNull()) {
          cweight += uC[*it - 1];
          cconvweight += uC[zC[*it - 1] - 1];
          cmultweight += uC[*it - 1] * uC[zC[*it - 1] - 1];
        }
      } else {
        idxhigh.push_back(*it);
        if (keyC == 1) {
          // "ken"
          noninv += cnum;
        } else if (keyC == 2) {
          // "aken"
          noninv += cmin;
        } else if (keyC == 3) {
          // "top"
          noninv += (((*it >= kC) && (zC[*it - 1] >= kC)) ? ctop : 0.0);
        } else if (keyC == 4) {
          // "add"
          noninv += cmultweight + cweight * uC[zC[*it - 1] - 1] + cconvweight * uC[*it - 1] + cnum * uC[*it - 1] * uC[zC[*it - 1] - 1];
        } else if (keyC == 5) {
          // "mult"
          noninv += cmultweight * uC[*it - 1] * uC[zC[*it - 1] - 1];
        }
      }
    }
    
    kendall_weight_quick(wrap(idxhigh), zC, nC, keyC, kC, uC);
    kendall_weight_quick(wrap(idxlow), zC, nC, keyC, kC, uC);
  }
  
  return noninv;
}
