// Copyright (c) 2025 genpca contributors
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// ---- helpers ---------------------------------------------------------------

static inline arma::uvec topk_indices_desc(const arma::vec& eval, const int k) {
  arma::uvec ord = arma::sort_index(eval, "descend");
  return ord.head(k);
}

// Try ARPACK top-k; fall back to full eigen if unavailable.
static bool eigs_topk(const arma::mat& M, const int k, arma::vec& eval, arma::mat& evec) {
  if (k <= 0 || k >= (int)M.n_rows) return false;
#ifdef ARMA_USE_ARPACK
  try {
    arma::eigs_sym(eval, evec, M, k, "la"); // largest algebraic
    return true;
  } catch (...) {
    return false;
  }
#else
  (void)M; (void)k; (void)eval; (void)evec;
  return false;
#endif
}

// Filter and sqrt eigenvalues -> singular values
static arma::vec sqrt_pos(const arma::vec& x, const double tol) {
  arma::vec y = x;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    y(i) = (y(i) > tol) ? std::sqrt(y(i)) : 0.0;
  }
  return y;
}

// ---- PRIMAL path (p <= n): needs Q and L_R (lower), returns scores/components ----
template <typename MatQ>
Rcpp::List gmd_primal_impl(const arma::mat& X,
                           const MatQ& Q,
                           const arma::mat& L_R,
                           const int k,
                           const double tol,
                           const bool topk) {
  const int p = (int)X.n_cols;
  const int k_use = std::min(k, p);
  // S = X' Q X
  arma::mat S = X.t() * (Q * X);

  // M = L_R^{-1} S L_R^{-T}
  arma::mat Linv = arma::inv(arma::trimatl(L_R));
  arma::mat M = Linv * S * Linv.t();

  arma::vec eval;
  arma::mat Z;
  bool ok = false;
  if (topk) ok = eigs_topk(M, k_use, eval, Z);
  if (!ok) {
    if (!arma::eig_sym(eval, Z, M)) Rcpp::stop("eig_sym failed (primal).");
    arma::uvec ord = topk_indices_desc(eval, k_use);
    Z = Z.cols(ord);
    eval = eval.elem(ord);
  }

  arma::vec d = sqrt_pos(eval, tol);

  // components C = R V = L_R^T Z
  arma::mat C = L_R.t() * Z;                // p x k

  // scores U = Q X C
  arma::mat U = (Q * X) * C;                // n x k

  return Rcpp::List::create(
      Rcpp::Named("u") = U,
      Rcpp::Named("v") = C,
      Rcpp::Named("d") = d
  );
}

// ---- DUAL path (n < p): needs L_Q (lower) and R, returns scores/components ----
template <typename MatR>
Rcpp::List gmd_dual_impl(const arma::mat& X,
                         const arma::mat& L_Q,
                         const MatR& R,
                         const int k,
                         const double tol,
                         const bool topk) {
  const int n = (int)X.n_rows;
  const int k_use = std::min(k, n);

  // B = L_Q^T X  (n x p)
  arma::mat B = L_Q.t() * X;

  // M = Q^{1/2} X R X^T Q^{1/2} = B R B^T  (n x n)
  arma::mat RBt = R * B.t();             // (p x n)
  arma::mat M   = B * RBt;               // (n x n)

  arma::vec eval;
  arma::mat Z;                           // Z = \tilde U (n x k)
  bool ok = false;
  if (topk) ok = eigs_topk(M, k_use, eval, Z);
  if (!ok) {
    if (!arma::eig_sym(eval, Z, M)) Rcpp::stop("eig_sym failed (dual).");
    arma::uvec ord = topk_indices_desc(eval, k_use);
    Z = Z.cols(ord);
    eval = eval.elem(ord);
  }
  arma::vec d = sqrt_pos(eval, tol);     // singular values

  // V = X^T (L_Q Z) / d   (p x k); Components C = R V
  arma::mat LQZ = L_Q * Z;               // (n x k)
  arma::mat Vtmp = X.t() * LQZ;          // (p x k)
  for (int i = 0; i < (int)d.n_elem; ++i) {
    if (d(i) > tol) Vtmp.col(i) /= d(i); else Vtmp.col(i).zeros();
  }
  arma::mat C = R * Vtmp;                // components = R V

  // Scores U = Q X C = (L_Q L_Q^T) X C = L_Q * (L_Q^T X) * C = L_Q * B * C
  arma::mat U = L_Q * (B * C);

  return Rcpp::List::create(
      Rcpp::Named("u") = U,
      Rcpp::Named("v") = C,
      Rcpp::Named("d") = d
  );
}

// ---- Exported entry points -------------------------------------------------

// Non-cached dense path: compute L_R / L_Q internally and dispatch primal/dual by size
template <typename MatQ, typename MatR>
Rcpp::List gmd_fast_auto(const arma::mat& X,
                         const MatQ& Q,
                         const MatR& R,
                         const int k,
                         const double tol,
                         const bool topk) {
  const int n = (int)X.n_rows;
  const int p = (int)X.n_cols;
  if (p <= n) {
    arma::mat Rdense(R);
    arma::mat L;
    if (!arma::chol(L, Rdense, "lower")) Rcpp::stop("Cholesky of R failed.");
    return gmd_primal_impl(X, Q, L, k, tol, topk);
  } else {
    arma::mat Qdense(Q);
    arma::mat L;
    if (!arma::chol(L, Qdense, "lower")) Rcpp::stop("Cholesky of Q failed.");
    return gmd_dual_impl(X, L, R, k, tol, topk);
  }
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dn(const arma::mat& X,
                           const arma::mat& Q,
                           const arma::mat& R,
                           const int k,
                           const double tol = 1e-8,
                           const bool topk = true) {
  return gmd_fast_auto(X, Q, R, k, tol, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_sp(const arma::mat& X,
                           const arma::sp_mat& Q,
                           const arma::sp_mat& R,
                           const int k,
                           const double tol = 1e-8,
                           const bool topk = true) {
  return gmd_fast_auto(X, Q, R, k, tol, topk);
}

// Cached primal/dual entry points (dense L factors provided by R)
// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_primal_dn(const arma::mat& X,
                                  const arma::mat& Q,
                                  const arma::mat& L_R,
                                  const int k,
                                  const double tol = 1e-8,
                                  const bool topk = true) {
  return gmd_primal_impl(X, Q, L_R, k, tol, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_primal_sp(const arma::mat& X,
                                  const arma::sp_mat& Q,
                                  const arma::mat& L_R,
                                  const int k,
                                  const double tol = 1e-8,
                                  const bool topk = true) {
  return gmd_primal_impl(X, Q, L_R, k, tol, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dual_dn(const arma::mat& X,
                                const arma::mat& L_Q,
                                const arma::mat& R,
                                const int k,
                                const double tol = 1e-8,
                                const bool topk = true) {
  return gmd_dual_impl(X, L_Q, R, k, tol, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dual_sp(const arma::mat& X,
                                const arma::mat& L_Q,
                                const arma::sp_mat& R,
                                const int k,
                                const double tol = 1e-8,
                                const bool topk = true) {
  return gmd_dual_impl(X, L_Q, R, k, tol, topk);
}