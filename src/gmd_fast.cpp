// Copyright (c) 2025 genpca contributors
#include <RcppArmadillo.h>
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <algorithm>
#include <cstring>
// [[Rcpp::depends(RcppArmadillo, RcppEigen, RSpectra)]]

// ---- helpers ---------------------------------------------------------------

static inline arma::uvec topk_indices_desc(const arma::vec& eval, const int k) {
  arma::uvec ord = arma::sort_index(eval, "descend");
  return ord.head(k);
}

static inline int choose_ncv(const int dim, const int k) {
  // Standard Lanczos heuristic: ncv > k and not too small.
  int ncv = std::max(2 * k + 1, 20);
  ncv = std::min(dim, ncv);
  if (ncv <= k) ncv = std::min(dim, k + 2);
  return ncv;
}

// Filter and sqrt eigenvalues -> singular values
static arma::vec sqrt_pos(const arma::vec& x, const double tol) {
  arma::vec y = x;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    y(i) = (y(i) > tol) ? std::sqrt(y(i)) : 0.0;
  }
  return y;
}

template <typename OpType>
static bool spectra_topk(OpType& op,
                         const int dim,
                         const int k,
                         const int maxit,
                         const double tol,
                         arma::vec& eval,
                         arma::mat& evec) {
  if (k <= 0 || k >= dim) return false;

  const int ncv = choose_ncv(dim, k);
  Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, OpType> eigs(&op, k, ncv);
  eigs.init();
  const int nconv = eigs.compute(maxit, tol, Spectra::LARGEST_ALGE);

  if (eigs.info() != Spectra::SUCCESSFUL || nconv != k) return false;

  const Eigen::VectorXd eval_e = eigs.eigenvalues();
  const Eigen::MatrixXd evec_e = eigs.eigenvectors();
  if (eval_e.size() == 0 || evec_e.cols() == 0) return false;

  eval.set_size(eval_e.size());
  std::memcpy(eval.memptr(), eval_e.data(), sizeof(double) * static_cast<size_t>(eval_e.size()));

  evec.set_size(evec_e.rows(), evec_e.cols());
  std::memcpy(evec.memptr(), evec_e.data(), sizeof(double) * static_cast<size_t>(evec_e.size()));
  return true;
}

// ---- matrix-free symmetric operators for Spectra ---------------------------

// Primal operator for M = L_R^{-1} (X' Q X) L_R^{-T}
template <typename MatQ>
class PrimalSymOp {
 public:
  PrimalSymOp(const arma::mat& X, const MatQ& Q, const arma::mat& L_R)
      : X_(X), Q_(Q), L_R_(L_R), dim_(static_cast<int>(X.n_cols)) {}

  int rows() const { return dim_; }
  int cols() const { return dim_; }

  void perform_op(const double* x_in, double* y_out) const {
    arma::vec x(const_cast<double*>(x_in), dim_, false, true);
    // z = L_R^{-T} x
    arma::vec z = arma::solve(arma::trimatu(L_R_.t()), x, arma::solve_opts::fast);
    // t = X' Q X z
    arma::vec t = X_ * z;
    t = Q_ * t;
    t = X_.t() * t;
    // y = L_R^{-1} t
    arma::vec y = arma::solve(arma::trimatl(L_R_), t, arma::solve_opts::fast);
    std::memcpy(y_out, y.memptr(), sizeof(double) * static_cast<size_t>(dim_));
  }

 private:
  const arma::mat& X_;
  const MatQ& Q_;
  const arma::mat& L_R_;
  const int dim_;
};

// Dual operator for M = L_Q^T X R X^T L_Q
template <typename MatR>
class DualSymOp {
 public:
  DualSymOp(const arma::mat& X, const arma::mat& L_Q, const MatR& R)
      : X_(X), L_Q_(L_Q), R_(R), dim_(static_cast<int>(X.n_rows)) {}

  int rows() const { return dim_; }
  int cols() const { return dim_; }

  void perform_op(const double* x_in, double* y_out) const {
    arma::vec x(const_cast<double*>(x_in), dim_, false, true);
    // y = L_Q^T X R X^T L_Q x
    arma::vec t = L_Q_ * x;
    t = X_.t() * t;
    t = R_ * t;
    t = X_ * t;
    arma::vec y = L_Q_.t() * t;
    std::memcpy(y_out, y.memptr(), sizeof(double) * static_cast<size_t>(dim_));
  }

 private:
  const arma::mat& X_;
  const arma::mat& L_Q_;
  const MatR& R_;
  const int dim_;
};

// ---- PRIMAL path (p <= n): needs Q and L_R (lower), returns scores/components ----
template <typename MatQ>
Rcpp::List gmd_primal_impl(const arma::mat& X,
                           const MatQ& Q,
                           const arma::mat& L_R,
                           const int k,
                           const double tol,
                           const int maxit,
                           const bool topk) {
  const int p = static_cast<int>(X.n_cols);
  const int k_use = std::min(k, p);

  arma::vec eval;
  arma::mat Z;
  bool ok = false;

  if (topk && k_use < p) {
    PrimalSymOp<MatQ> op(X, Q, L_R);
    ok = spectra_topk(op, p, k_use, maxit, tol, eval, Z);
  }

  if (!ok) {
    // Fallback: full matrix construction + dense eigendecomposition.
    arma::mat S = X.t() * (Q * X);
    arma::mat Linv = arma::inv(arma::trimatl(L_R));
    arma::mat M = Linv * S * Linv.t();
    if (!arma::eig_sym(eval, Z, M)) Rcpp::stop("eig_sym failed (primal fallback).");
    arma::uvec ord = topk_indices_desc(eval, k_use);
    Z = Z.cols(ord);
    eval = eval.elem(ord);
  }

  if (eval.n_elem == 0 || Z.n_cols == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec()
    );
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
                         const int maxit,
                         const bool topk) {
  const int n = static_cast<int>(X.n_rows);
  const int k_use = std::min(k, n);

  arma::vec eval;
  arma::mat Z;                           // Z = \tilde U (n x k)
  bool ok = false;

  if (topk && k_use < n) {
    DualSymOp<MatR> op(X, L_Q, R);
    ok = spectra_topk(op, n, k_use, maxit, tol, eval, Z);
  }

  if (!ok) {
    // Fallback: full matrix construction + dense eigendecomposition.
    arma::mat B = L_Q.t() * X;
    arma::mat RBt = R * B.t();
    arma::mat M = B * RBt;
    if (!arma::eig_sym(eval, Z, M)) Rcpp::stop("eig_sym failed (dual fallback).");
    arma::uvec ord = topk_indices_desc(eval, k_use);
    Z = Z.cols(ord);
    eval = eval.elem(ord);
  }

  if (eval.n_elem == 0 || Z.n_cols == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec()
    );
  }

  arma::vec d = sqrt_pos(eval, tol);     // singular values

  // V = X^T (L_Q Z) / d   (p x k); Components C = R V
  arma::mat LQZ = L_Q * Z;               // (n x k)
  arma::mat Vtmp = X.t() * LQZ;          // (p x k)
  for (int i = 0; i < static_cast<int>(d.n_elem); ++i) {
    if (d(i) > tol) Vtmp.col(i) /= d(i); else Vtmp.col(i).zeros();
  }
  arma::mat C = R * Vtmp;                // components = R V

  // Scores U = Q X C = (L_Q L_Q^T) X C = L_Q * (L_Q^T X) * C
  arma::mat U = L_Q * ((L_Q.t() * X) * C);

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
                         const int maxit,
                         const bool topk) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  if (p <= n) {
    arma::mat Rdense(R);
    arma::mat L;
    if (!arma::chol(L, Rdense, "lower")) Rcpp::stop("Cholesky of R failed.");
    return gmd_primal_impl(X, Q, L, k, tol, maxit, topk);
  } else {
    arma::mat Qdense(Q);
    arma::mat L;
    if (!arma::chol(L, Qdense, "lower")) Rcpp::stop("Cholesky of Q failed.");
    return gmd_dual_impl(X, L, R, k, tol, maxit, topk);
  }
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dn(const arma::mat& X,
                           const arma::mat& Q,
                           const arma::mat& R,
                           const int k,
                           const double tol = 1e-8,
                           const int maxit = 1000,
                           const bool topk = true) {
  return gmd_fast_auto(X, Q, R, k, tol, maxit, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_sp(const arma::mat& X,
                           const arma::sp_mat& Q,
                           const arma::sp_mat& R,
                           const int k,
                           const double tol = 1e-8,
                           const int maxit = 1000,
                           const bool topk = true) {
  return gmd_fast_auto(X, Q, R, k, tol, maxit, topk);
}

// Cached primal/dual entry points (dense L factors provided by R)
// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_primal_dn(const arma::mat& X,
                                  const arma::mat& Q,
                                  const arma::mat& L_R,
                                  const int k,
                                  const double tol = 1e-8,
                                  const int maxit = 1000,
                                  const bool topk = true) {
  return gmd_primal_impl(X, Q, L_R, k, tol, maxit, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_primal_sp(const arma::mat& X,
                                  const arma::sp_mat& Q,
                                  const arma::mat& L_R,
                                  const int k,
                                  const double tol = 1e-8,
                                  const int maxit = 1000,
                                  const bool topk = true) {
  return gmd_primal_impl(X, Q, L_R, k, tol, maxit, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dual_dn(const arma::mat& X,
                                const arma::mat& L_Q,
                                const arma::mat& R,
                                const int k,
                                const double tol = 1e-8,
                                const int maxit = 1000,
                                const bool topk = true) {
  return gmd_dual_impl(X, L_Q, R, k, tol, maxit, topk);
}

// [[Rcpp::export]]
Rcpp::List gmd_fast_cpp_dual_sp(const arma::mat& X,
                                const arma::mat& L_Q,
                                const arma::sp_mat& R,
                                const int k,
                                const double tol = 1e-8,
                                const int maxit = 1000,
                                const bool topk = true) {
  return gmd_dual_impl(X, L_Q, R, k, tol, maxit, topk);
}
