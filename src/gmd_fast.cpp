// Copyright (c) 2025 genpca contributors
#include <RcppArmadillo.h>
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <random>
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

static inline arma::mat random_sign_matrix(const arma::uword nrow,
                                           const arma::uword ncol,
                                           const unsigned int seed) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> d01(0, 1);
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = d01(gen) ? 1.0 : -1.0;
    }
  }
  return out;
}

template <typename MatM>
static arma::mat metric_orthonormalize_cpp(const arma::mat& A,
                                           const MatM& M,
                                           const double jitter,
                                           const double tol) {
  if (A.n_cols == 0) {
    return arma::mat(A.n_rows, 0, arma::fill::zeros);
  }

  arma::mat MA = M * A;
  arma::mat G = A.t() * MA;
  G = 0.5 * (G + G.t());

  double jitter_now = std::max(0.0, jitter);
  for (int tries = 0; tries < 7; ++tries) {
    arma::mat Greg = G;
    if (jitter_now > 0.0) {
      Greg.diag() += jitter_now;
    }
    arma::mat C;
    if (arma::chol(C, Greg, "upper")) {
      // A * inv(C), implemented as triangular solve.
      arma::mat Xt = arma::solve(arma::trimatl(C.t()), A.t(), arma::solve_opts::fast);
      return Xt.t();
    }
    jitter_now = (jitter_now > 0.0) ? jitter_now * 10.0 : 1e-10;
  }

  arma::vec eval;
  arma::mat evec;
  if (!arma::eig_sym(eval, evec, G)) {
    return arma::mat(A.n_rows, 0, arma::fill::zeros);
  }
  arma::uvec keep = arma::find(eval > tol);
  if (keep.n_elem == 0) {
    return arma::mat(A.n_rows, 0, arma::fill::zeros);
  }

  arma::vec invsqrt = 1.0 / arma::sqrt(eval.elem(keep));
  arma::mat B = evec.cols(keep) * arma::diagmat(invsqrt);
  return A * B;
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

// Primal operator for M = L_R^T (X' Q X) L_R
template <typename MatQ>
class PrimalSymOp {
 public:
  PrimalSymOp(const arma::mat& X, const MatQ& Q, const arma::mat& L_R)
      : X_(X), Q_(Q), L_R_(L_R), dim_(static_cast<int>(X.n_cols)) {}

  int rows() const { return dim_; }
  int cols() const { return dim_; }

  void perform_op(const double* x_in, double* y_out) const {
    arma::vec x(const_cast<double*>(x_in), dim_, false, true);
    // z = L_R x
    arma::vec z = L_R_ * x;
    // t = X' Q X z
    arma::vec t = X_ * z;
    t = Q_ * t;
    t = X_.t() * t;
    // y = L_R^T t
    arma::vec y = L_R_.t() * t;
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

// ---- PRIMAL path (p <= n): needs Q and L_R (lower), returns scores/components + ou/ov ----
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
    arma::mat M = L_R.t() * S * L_R;
    if (!arma::eig_sym(eval, Z, M)) Rcpp::stop("eig_sym failed (primal fallback).");
    arma::uvec ord = topk_indices_desc(eval, k_use);
    Z = Z.cols(ord);
    eval = eval.elem(ord);
  }

  if (eval.n_elem == 0 || Z.n_cols == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("ou") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("ov") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec()
    );
  }

  arma::vec d = sqrt_pos(eval, tol);

  // ov = R^{-1/2} Z where Z are eigenvectors of L_R^T (X'QX) L_R.
  arma::mat ov = arma::solve(arma::trimatu(L_R.t()), Z, arma::solve_opts::fast);

  // Components C = R ov = L_R Z
  arma::mat C = L_R * Z;                    // p x k

  // ou from X R ov = ou diag(d)
  arma::mat ou = X * C;                     // n x k
  for (int i = 0; i < static_cast<int>(d.n_elem); ++i) {
    if (d(i) > tol) ou.col(i) /= d(i); else ou.col(i).zeros();
  }

  // scores U = Q X C
  arma::mat U = (Q * X) * C;                // n x k

  return Rcpp::List::create(
      Rcpp::Named("u") = U,
      Rcpp::Named("v") = C,
      Rcpp::Named("ou") = ou,
      Rcpp::Named("ov") = ov,
      Rcpp::Named("d") = d
  );
}

// ---- DUAL path (n < p): needs L_Q (lower) and R, returns scores/components + ou/ov ----
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
      Rcpp::Named("ou") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("ov") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec()
    );
  }

  arma::vec d = sqrt_pos(eval, tol);     // singular values

  // ou = Q^{-1/2} Z where Z are eigenvectors of L_Q^T (X R X^T) L_Q.
  arma::mat ou = arma::solve(arma::trimatu(L_Q.t()), Z, arma::solve_opts::fast);

  // ov = X^T Q ou / d = X^T (L_Q Z) / d
  arma::mat ov = X.t() * (L_Q * Z);      // (p x k)
  for (int i = 0; i < static_cast<int>(d.n_elem); ++i) {
    if (d(i) > tol) ov.col(i) /= d(i); else ov.col(i).zeros();
  }

  // Components C = R ov
  arma::mat C = R * ov;

  // Scores U = Q X C = (L_Q L_Q^T) X C = L_Q * (L_Q^T X) * C
  arma::mat U = L_Q * ((L_Q.t() * X) * C);

  return Rcpp::List::create(
      Rcpp::Named("u") = U,
      Rcpp::Named("v") = C,
      Rcpp::Named("ou") = ou,
      Rcpp::Named("ov") = ov,
      Rcpp::Named("d") = d
  );
}

template <typename MatQ, typename MatR>
static void randomized_polish_cpp(const arma::mat& X,
                                  const MatQ& Q,
                                  const MatR& R,
                                  arma::mat& U,
                                  arma::mat& V,
                                  arma::vec& d,
                                  const int n_polish,
                                  const double jitter,
                                  const double tol,
                                  const double polish_tol) {
  if (n_polish <= 0 || U.n_cols == 0) return;

  arma::vec d_prev;
  for (int it = 0; it < n_polish; ++it) {
    arma::mat Y = X * (R * V);
    U = metric_orthonormalize_cpp(Y, Q, jitter, tol);
    arma::mat Z = X.t() * (Q * U);
    V = metric_orthonormalize_cpp(Z, R, jitter, tol);
    arma::mat RV = R * V;
    // T = U^T Q X R V; reuse Z = X^T Q U to avoid another pass through X.
    arma::mat T = Z.t() * RV;

    arma::mat P, Vrot;
    arma::vec s;
    if (!arma::svd(P, s, Vrot, T)) {
      Rcpp::stop("SVD failed in randomized polish.");
    }
    U = U * P;
    V = V * Vrot;
    d = s;

    if (polish_tol > 0.0 && d_prev.n_elem == d.n_elem && d.n_elem > 0) {
      arma::vec denom = arma::max(arma::abs(d_prev), arma::vec(d_prev.n_elem, arma::fill::value(1e-12)));
      arma::vec rel = arma::abs(d - d_prev) / denom;
      double max_rel = rel.max();
      if (std::isfinite(max_rel) && max_rel < polish_tol) {
        break;
      }
    }
    d_prev = d;
  }
}

template <typename MatQ, typename MatR>
Rcpp::List gmd_randomized_impl(const arma::mat& X,
                               const MatQ& Q,
                               const MatR& R,
                               const int k,
                               const int oversample,
                               const int n_power,
                               const int n_polish,
                               const double jitter,
                               const double tol,
                               const double polish_tol,
                               const unsigned int seed) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  const int k_use = std::max(0, std::min(k, std::min(n, p)));

  if (k_use == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }

  const int ell = std::max(1, std::min(std::min(n, p), k_use + std::max(0, oversample)));
  arma::mat Omega = random_sign_matrix(static_cast<arma::uword>(p),
                                       static_cast<arma::uword>(ell),
                                       seed);

  arma::mat Y = X * (R * Omega);
  for (int it = 0; it < n_power; ++it) {
    arma::mat Utmp = metric_orthonormalize_cpp(Y, Q, jitter, tol);
    if (Utmp.n_cols == 0) break;
    arma::mat Z = X.t() * (Q * Utmp);
    Y = X * (R * Z);
  }

  arma::mat U0 = metric_orthonormalize_cpp(Y, Q, jitter, tol);
  if (U0.n_cols == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }

  arma::mat B = X.t() * (Q * U0);
  arma::mat RB = R * B;
  arma::mat G = B.t() * RB;
  G = 0.5 * (G + G.t());

  arma::vec eval;
  arma::mat Sfull;
  if (!arma::eig_sym(eval, Sfull, G)) {
    Rcpp::stop("eig_sym failed in randomized solver.");
  }
  arma::uvec ord = arma::sort_index(eval, "descend");
  const arma::uword kk = std::min<arma::uword>(static_cast<arma::uword>(k_use), ord.n_elem);
  if (kk == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }
  arma::uvec take = ord.head(kk);

  arma::vec lam = eval.elem(take);
  for (arma::uword i = 0; i < lam.n_elem; ++i) {
    if (lam(i) < 0.0) lam(i) = 0.0;
  }
  arma::vec d = arma::sqrt(lam);
  arma::mat S = Sfull.cols(take);

  arma::mat U = U0 * S;
  arma::mat V = B * S;
  for (arma::uword i = 0; i < d.n_elem; ++i) {
    if (d(i) > tol) {
      V.col(i) /= d(i);
    } else {
      V.col(i).zeros();
    }
  }

  randomized_polish_cpp(X, Q, R, U, V, d, n_polish, jitter, tol, polish_tol);

  arma::uvec keep = arma::find(d > tol);
  if (keep.n_elem == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }

  U = U.cols(keep);
  V = V.cols(keep);
  d = d.elem(keep);

  return Rcpp::List::create(
    Rcpp::Named("u") = U,
    Rcpp::Named("v") = V,
    Rcpp::Named("d") = d,
    Rcpp::Named("k") = static_cast<int>(d.n_elem)
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

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_dn(const arma::mat& X,
                                 const arma::mat& Q,
                                 const arma::mat& R,
                                 const int k,
                                 const int oversample = 20,
                                 const int n_power = 1,
                                 const int n_polish = 0,
                                 const double jitter = 1e-10,
                                 const double tol = 1e-9,
                                 const double polish_tol = 0.0,
                                 const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_sp(const arma::mat& X,
                                 const arma::sp_mat& Q,
                                 const arma::sp_mat& R,
                                 const int k,
                                 const int oversample = 20,
                                 const int n_power = 1,
                                 const int n_polish = 0,
                                 const double jitter = 1e-10,
                                 const double tol = 1e-9,
                                 const double polish_tol = 0.0,
                                 const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_qsp_rdn(const arma::mat& X,
                                      const arma::sp_mat& Q,
                                      const arma::mat& R,
                                      const int k,
                                      const int oversample = 20,
                                      const int n_power = 1,
                                      const int n_polish = 0,
                                      const double jitter = 1e-10,
                                      const double tol = 1e-9,
                                      const double polish_tol = 0.0,
                                      const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_qdn_rsp(const arma::mat& X,
                                      const arma::mat& Q,
                                      const arma::sp_mat& R,
                                      const int k,
                                      const int oversample = 20,
                                      const int n_power = 1,
                                      const int n_polish = 0,
                                      const double jitter = 1e-10,
                                      const double tol = 1e-9,
                                      const double polish_tol = 0.0,
                                      const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}
