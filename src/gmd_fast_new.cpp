// Copyright (c) 2025 genpca contributors
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Solve generalized PCA via eigen on (X^T Q X, R), return top-k.
// Templated on Q and R to support arma::mat and arma::sp_mat.
template <typename MatQ, typename MatR>
Rcpp::List gmd_fast_impl(const arma::mat& X,
                         const MatQ& Q,
                         const MatR& R,
                         const int k,
                         const double tol) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int k_use = std::min(k, std::min(n, p));
  
  // Convert Q and R to dense for uniform processing
  arma::mat Qd = arma::symmatu(arma::mat(Q));
  arma::mat Rd = arma::symmatu(arma::mat(R));
  
  // A = X^T Q X (symmetric p x p)
  arma::mat A = X.t() * (Qd * X);
  
  // Generalized eigenproblem: A v = lambda R v
  // Transform to standard eigenproblem via Cholesky of R
  arma::vec eval;
  arma::mat V;
  
  // Cholesky decomposition of R: R = L L^T
  arma::mat L;
  bool chol_ok = arma::chol(L, Rd, "lower");
  if (!chol_ok) {
    Rcpp::stop("Cholesky decomposition of R failed (R not SPD?).");
  }
  
  // Transform to standard eigenproblem: L^{-1} A L^{-T} z = lambda z
  arma::mat Linv = arma::inv(arma::trimatl(L));
  arma::mat B = Linv * A * Linv.t();
  
  // Solve standard eigenproblem
  bool eig_ok = arma::eig_sym(eval, V, B);
  if (!eig_ok) {
    Rcpp::stop("Eigenvalue decomposition failed.");
  }
  
  // Back-transform eigenvectors: v = L^{-T} z
  V = Linv.t() * V;
  
  // Take top-k eigenpairs (largest eigenvalues)
  arma::uvec ord = arma::sort_index(eval, "descend");
  if (ord.n_elem > (size_t)k_use) {
    ord = ord.head(k_use);
  }
  arma::vec lam = eval.elem(ord);
  arma::mat Vk = V.cols(ord);
  
  // Singular values are sqrt(lambda) (clamp negatives due to numeric noise)
  lam.transform([tol](double v) { return (v > tol) ? std::sqrt(v) : 0.0; });
  
  // U = X V / d
  arma::mat U = X * Vk;
  for (int i = 0; i < k_use; ++i) {
    double d = lam(i);
    if (d > tol) {
      U.col(i) /= d;
    } else {
      U.col(i).zeros();
    }
  }
  
  // Q-orthonormalize U (modified Gram-Schmidt under Q)
  for (int j = 0; j < k_use; ++j) {
    for (int i = 0; i < j; ++i) {
      double alpha = arma::as_scalar(U.col(i).t() * (Qd * U.col(j)));
      U.col(j) -= alpha * U.col(i);
    }
    double nrm2 = arma::as_scalar(U.col(j).t() * (Qd * U.col(j)));
    double nrm = std::sqrt(std::max(tol, nrm2));
    if (nrm > tol) {
      U.col(j) /= nrm;
    }
  }
  
  // R-orthonormalize V (optional but helps comparisons)
  for (int j = 0; j < k_use; ++j) {
    for (int i = 0; i < j; ++i) {
      double alpha = arma::as_scalar(Vk.col(i).t() * (Rd * Vk.col(j)));
      Vk.col(j) -= alpha * Vk.col(i);
    }
    double nrm2 = arma::as_scalar(Vk.col(j).t() * (Rd * Vk.col(j)));
    double nrm = std::sqrt(std::max(tol, nrm2));
    if (nrm > tol) {
      Vk.col(j) /= nrm;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("u") = U,
    Rcpp::Named("v") = Vk,
    Rcpp::Named("d") = lam,
    Rcpp::Named("k") = k_use
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dn(const arma::mat& X,
                           const arma::mat& Q,
                           const arma::mat& R,
                           const int k,
                           const double tol = 1e-8) {
  return gmd_fast_impl(X, Q, R, k, tol);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_sp(const arma::mat& X,
                           const arma::sp_mat& Q,
                           const arma::sp_mat& R,
                           const int k,
                           const double tol = 1e-8) {
  return gmd_fast_impl(X, Q, R, k, tol);
}