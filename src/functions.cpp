// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

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
  double SCI = (N/W) * (cv/v) * (0.5);

  return SCI;
}


// [[Rcpp::export]]
arma::mat spatialCrossCorMatrix_C(arma::mat sigMat, arma::mat weight, bool display_progress=true) {
  int N = sigMat.n_rows;

  arma::mat scor;
  scor.set_size(N, N);

  int i;
  int j;
  arma::vec x;
  arma::vec y;
  double SCI;

  Progress p(N*N, display_progress);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (Progress::check_abort() ) {
        return scor;
      }
      p.increment();
      x = conv_to<arma::vec>::from(sigMat.row(i)); // convert from rowvec to vec
      y = conv_to<arma::vec>::from(sigMat.row(j));
      SCI = spatialCrossCor_C(x, y, weight);
      scor.at(i,j) = SCI;
    }
  }

  return scor;
}


// [[Rcpp::export]]
arma::vec moranTest_C(arma::vec x, arma::mat weight) {
  double N = weight.n_rows;

  // first moment
  double ei = -1/(N - 1);

  // unitization
  weight = unitize_C(weight);

  // Moran's I
  double W = sum(sum(weight));
  arma::vec z = x - mean(x);
  double cv = sum(sum(weight % (z * z.t()))); // weight * (z %o% z)
  NumericVector zbar = wrap(z); // convert to numeric vector to use power
  double v = sum(pow(zbar, 2));
  double obs = (N/W) * (cv/v);

  // second moment
  double Wsq = pow(W, 2);
  double Nsq = pow(N, 2);
  NumericMatrix wbar = wrap(weight + weight.t());
  arma::mat wbarbar = pow(wbar, 2);
  double S1 = 0.5 * sum(sum(wbarbar));
  arma::vec rs = conv_to<arma::vec>::from(sum(weight, 1));
  arma::vec cs = conv_to<arma::vec>::from(sum(weight, 0));
  arma::vec sg =  rs + cs;
  NumericVector sgbar = wrap(sg);
  arma::vec sbarbar = pow(sgbar, 2);
  double S2 = sum(sbarbar); //sg^2
  arma::vec zbarbar = pow(zbar, 4);
  double S3 = (sum(zbarbar)/N)/pow(v/N, 2);
  double S4 = (Nsq - 3*N + 3)*S1 - N*S2 + 3*Wsq;
  double S5 = (Nsq - N)*S1 - 2*N*S2 + 6*Wsq;
  double ei2 = (N*S4 - S3*S5)/((N - 1)*(N - 2)*(N - 3) * Wsq);

  // standard deviation
  double sdi = sqrt(ei2 - pow(ei, 2));

  // return results as vector
  arma::vec results;
  results.set_size(3);
  results.at(0) = obs; // Moran's I
  results.at(1) = ei; // Expected
  results.at(2) = sdi; // SD

  return results;
}


// [[Rcpp::export]]
arma::mat getSpatialPatterns_C(arma::mat mat, arma::mat adj, bool display_progress=true) {
  double N = mat.n_rows;

  arma::mat results;
  results.set_size(N, 3);

  int i;
  arma::vec value;
  arma::vec Ir;

  Progress p(N, display_progress);
  for (i = 0; i < N; i++) {
    if (Progress::check_abort() ) {
      return results;
    }
    p.increment();
    value = conv_to<arma::vec>::from(mat.row(i));
    Ir = moranTest_C(value, adj);
    results.row(i) = Ir.t();
  }

  return results;
}


