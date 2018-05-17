// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat unitize_C(arma::mat weight) {

  int N = weight.n_cols;

  // scale weights
  int i;
  int j;
  arma::vec rs = sum(weight, 1);
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

  return weight;
}


// [[Rcpp::export]]
double spatialCrossCor_C(arma::vec x,
                            arma::vec y,
                            arma::mat weight) {
  int N = weight.n_cols;
  weight = unitize_C(weight);

  // compute spatial cross correlation
  double W = sum(sum(weight));
  arma::vec dx = x - mean(x);
  arma::vec dy = y - mean(y);

  arma::mat cv1 = dx * dy.t(); // outer product
  arma::mat cv2 = dy * dx.t();

  double cv = sum(sum(weight % ( cv1 + cv2 ))); // element-wise product

  arma::vec v1 = pow(dx, 2);
  arma::vec v2 = pow(dy, 2);
  double v = sqrt(sum(v1) * sum(v2));
  double SCI = (N/W) * (cv/v);

  return SCI;
}


// [[Rcpp::export]]
arma::mat spatialCrossCorMatrix_C(arma::mat sigMat, arma::mat weight) {
  int N = sigMat.n_cols;
  int M = sigMat.n_rows;

  arma::mat scor;
  scor.set_size(N, M);

  int i;
  int j;
  arma::vec x;
  arma::vec y;
  double SCI;
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      x = conv_to<arma::vec>::from(sigMat.row(i)); // convert from rowvec to vec
      y = conv_to<arma::vec>::from(sigMat.row(j));
      SCI = spatialCrossCor_C(x, y, weight);
      scor.at(i,j) = SCI;
    }
  }

  return scor;
}


/*** R
w = matrix(c(1,1,0,1), 2,2)
#unitize_C(w)
#weight = w
#rs <- rowSums(weight)
#rs[rs == 0] <- 1
#weight <- weight/rs
#print(weight)

x = c(1,2)
y = c(1,3)
#spatialCrossCor_C(x,y,w)

sigMat = rbind(x,y)
spatialCrossCorMatrix_C(sigMat, w)
*/
