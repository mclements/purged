#include <Rcpp.h>

using namespace Rcpp;

class Rpexp {
public: 
  Rpexp() {} // blank default constructor
  Rpexp(double *hin, double *tin, int nin) : n(nin) {
    int i;
    H.resize(n);
    t.resize(n);
    h.resize(n);
    H[0]=0.0; h[0]=hin[0]; t[0]=tin[0];
    if (n>1) {
      for(i=1;i<n;i++) {
	h[i]=hin[i]; t[i]=tin[i];
	H[i] = H[i-1]+(t[i]-t[i-1])*h[i-1];
      }
    }
  }
  double rand(double from = 0.0) {
    double v = 0.0, H0 = 0.0, tstar = 0.0;
    int i = 0, i0 = 0;
    if (from > 0.0) {
      i0 = (from >= t[n-1]) ? (n-1) : int(std::lower_bound(t.begin(), t.end(), from) - t.begin())-1;
      H0 = H[i0] + (from - t[i0])*h[i0];
    }
    v = R::rexp(1.0) + H0;
    i = (v >= H[n-1]) ? (n-1) : int(std::lower_bound(H.begin(), H.end(), v) - H.begin())-1;
    tstar = t[i]+(v-H[i])/h[i];
    return tstar;
  }
 private:
  std::vector<double> H, h, t;
  int n;
};


// [[Rcpp::export]]
NumericVector 
rpexp(const int n, NumericVector t, NumericVector rate, NumericVector t0) {
  Rpexp random(&rate[0], &t[0], t.size());
  NumericVector out(n);
  for (int i=0; i<n; ++i) {
    out[i] = random.rand(t0[i]);
  }
  return out;
}

// R: .Call("purged_rpexp",10L,c(0,1),c(0.1,0.2),c(0,0),package="purged")
