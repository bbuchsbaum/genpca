#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Generalized power-iteration deflation (Allen, Grosenick & Taylor, 2014):
//   uhat = Xhat R v;  u = uhat / sqrt(u' Q u)
//   vhat = Xhat' Q u; v = vhat / sqrt(v' R v)
//   d_i  = u' Q X R v
// Degenerate directions (zero norm in the metric) or near-zero singular
// values stop extraction early and results are trimmed, matching the R
// implementation gmd_deflationR().
//
// The residual is applied implicitly.  After j components have been extracted:
//   X_j  w = X w  - U_j (d_j * (V_j' w))
//   X_j' z = X' z - V_j (d_j * (U_j' z))
// This keeps sparse X sparse and avoids a full dense residual copy.

// large dimension should be in ROWS (if columns is X, then pass (t(X), R,))

template <typename MatX>
List gmd_deflation_impl(const MatX &X, const arma::sp_mat &Q, const arma::sp_mat &R,
                        int k, double thr=1e-7, int maxit=500, bool verbose=false) {
  if (maxit < 1) {
    stop("maxit must be >= 1.");
  }
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat ugmd(n, k, fill::zeros);
  arma::mat vgmd(p, k, fill::zeros);
  arma::vec dgmd(k, fill::zeros);
  arma::vec propv(k, fill::zeros);
  arma::vec cumv(k, fill::zeros);

  arma::vec u = randn(n);
  arma::vec v = randn(p);

  // qrnorm = tr(X' Q X R) = accu((Q X) % (X R)); avoids any n x n / p x p
  // dense intermediate (trace is cyclic, so one expression covers both shapes).
  double qrnorm = arma::accu((Q * X) % (X * R));
  if (!std::isfinite(qrnorm) || qrnorm < 1e-12) {
    warning("Total generalized variance is near zero; explained variance may be unstable.");
    qrnorm = 1.0;
  }

  int k_found = 0;

  auto residual_mv = [&](const arma::vec& w, const int n_prev) -> arma::vec {
    arma::vec y = X * w;
    if (n_prev > 0) {
      const arma::uword nprev = static_cast<arma::uword>(n_prev);
      const arma::uword last = nprev - 1;
      arma::vec coeff = dgmd.head(nprev) % (vgmd.cols(0, last).t() * w);
      y -= ugmd.cols(0, last) * coeff;
    }
    return y;
  };

  auto residual_t_mv = [&](const arma::vec& z, const int n_prev) -> arma::vec {
    arma::vec y = X.t() * z;
    if (n_prev > 0) {
      const arma::uword nprev = static_cast<arma::uword>(n_prev);
      const arma::uword last = nprev - 1;
      arma::vec coeff = dgmd.head(nprev) % (ugmd.cols(0, last).t() * z);
      y -= vgmd.cols(0, last) * coeff;
    }
    return y;
  };

  for (int i=0; i<k; i++) {
    double err = 1;
    int iter = 0;
    bool degenerate = false;
    while (err > thr && iter < maxit) {
      iter += 1;
      arma::vec oldu = u;
      arma::vec oldv = v;

      arma::vec uhat = residual_mv(R * v, k_found);
      double u_norm = std::sqrt(arma::as_scalar(uhat.t() * (Q * uhat)));
      if (!std::isfinite(u_norm) || u_norm <= 0.0) { degenerate = true; break; }
      u = uhat / u_norm;

      arma::vec vhat = residual_t_mv(Q * u, k_found);
      double v_norm = std::sqrt(arma::as_scalar(vhat.t() * (R * vhat)));
      if (!std::isfinite(v_norm) || v_norm <= 0.0) { degenerate = true; break; }
      v = vhat / v_norm;

      err = arma::accu(arma::square(oldu - u)) + arma::accu(arma::square(oldv - v));
    }

    if (degenerate) {
      if (verbose) {
        Rcout << "Degenerate direction at component " << (i + 1) << "; stopping deflation." << std::endl;
      }
      warning("Deflation stopped early at component %d (degenerate direction).", i + 1);
      break;
    }

    if (iter >= maxit && err > thr) {
      if (verbose) {
        Rcout << "Power iteration reached maxit for component " << (i + 1) << "." << std::endl;
      }
      warning("Power iteration reached maxit at component %d (maxit=%d).", i + 1, maxit);
    }

    double d_i = arma::as_scalar((Q * u).t() * residual_mv(R * v, k_found));
    if (!std::isfinite(d_i) || std::fabs(d_i) < thr) {
      warning("Deflation stopped early at component %d (singular value near zero).", i + 1);
      break;
    }

    dgmd(k_found) = d_i;
    ugmd.col(k_found) = u;
    vgmd.col(k_found) = v;
    propv(k_found) = d_i * d_i / qrnorm;
    cumv(k_found) = (k_found > 0) ? (cumv(k_found - 1) + propv(k_found)) : propv(k_found);
    k_found += 1;
  }

  return List::create(
    Named("d") = dgmd.head(k_found),
    Named("v") = vgmd.head_cols(k_found),
    Named("u") = ugmd.head_cols(k_found),
    Named("k") = k_found,
    Named("cumv") = cumv.head(k_found),
    Named("propv") = propv.head(k_found)
  );
}

//[[Rcpp::export]]
List gmd_deflation_cpp(const arma::mat &X, const arma::sp_mat &Q, const arma::sp_mat &R,
                       int k, double thr=1e-7, int maxit=500, bool verbose=false) {
  return gmd_deflation_impl(X, Q, R, k, thr, maxit, verbose);
}

//[[Rcpp::export]]
List gmd_deflation_cpp_sp(const arma::sp_mat &X, const arma::sp_mat &Q, const arma::sp_mat &R,
                          int k, double thr=1e-7, int maxit=500, bool verbose=false) {
  return gmd_deflation_impl(X, Q, R, k, thr, maxit, verbose);
}
