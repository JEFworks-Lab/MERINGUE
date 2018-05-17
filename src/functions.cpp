// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double spatialCrossCor_C(arma::vec x,
                            arma::vec y,
                            arma::mat weight) {

  double N;
  double W;
  arma::vec dx;
  arma::vec dy;
  arma::vec rs;

  N = weight.n_cols;

  // scale weights
  int i;
  int j;
  rs = sum(weight, 1);
  for (i = 0; i < N; i++) {
    if(rs.at(i) == 0) {
      rs.at(i) = 1;
    }
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      weight.at(i,j) = weight.at(i,j) / rs.at(i);
    }
  }

  // compute spatial cross correlation
  W = sum(rs);
  dx = x - mean(x);
  dy = y - mean(y);

  arma::mat cv1;
  arma::mat cv2;
  cv1 = dx * dy.t();
  cv2 = dy * dx.t();

  double cv;
  cv = sum(sum(weight * ( cv1 + cv2 )));
  arma::vec v1 = pow(dx, 2);
  arma::vec v2 = pow(dy, 2);
  double v = sqrt(sum(v1) * sum(v2));
  double SCI = (N/W) * (cv/v);

  return SCI;
}

/*** R
spatialCrossCor_C(c(1,2), c(1,3), matrix(c(1,1,0,1), 2,2))
*/
