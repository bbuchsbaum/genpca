// [[Rcpp::depends(RcppArmadillo, RSpectra)]]
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD 1
#endif

#include <RcppArmadillo.h>
#include <RcppEigen.h>   

#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

using namespace Rcpp;
using namespace arma;

// ──────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────

// diagonal‑only SPD square root + inverse
// Checks if matrix is diagonal before proceeding.
inline void compute_diag_sqrt_inv(const sp_mat& D,
                                  sp_mat& Dhalf,
                                  sp_mat& Dhalf_inv,
                                  double tol = 1e-12)
{
    if (!D.is_diagmat()) {
        Rcpp::stop("compute_diag_sqrt_inv called on non-diagonal matrix.");
    }
    vec d_diag = vec(D.diag()); // Explicitly construct vec from spdiagview
    uvec idx = find(d_diag > tol);
    if (idx.n_elem == 0) {
        Rcpp::warning("No diagonal elements > tol found in compute_diag_sqrt_inv.");
        Dhalf = sp_mat(D.n_rows, D.n_cols);
        Dhalf_inv = sp_mat(D.n_rows, D.n_cols);
        return;
    }
    vec s = sqrt(d_diag(idx));

    Dhalf      = sp_mat(D.n_rows, D.n_cols);
    Dhalf_inv  = sp_mat(D.n_rows, D.n_cols);

    // Useumat constructor for sparse diagonal matrix
    umat locations(2, idx.n_elem);
    locations.row(0) = idx.t(); 
    locations.row(1) = idx.t();

    Dhalf     = sp_mat(locations, s, D.n_rows, D.n_cols);
    Dhalf_inv = sp_mat(locations, 1.0 / s, D.n_rows, D.n_cols);
}

// dense or sparse SPD square root + inverse (eigen route)
// Added size check.
inline void compute_spd_sqrt_inv(const sp_mat& A,
                                 sp_mat& Ahalf,
                                 sp_mat& Ahalf_inv,
                                 double tol = 1e-12,
                                 unsigned int max_dense_dim = 4000)
{
    if (A.n_rows > max_dense_dim) {
        Rcpp::stop("Constraint matrix too large (%d x %d) for dense eigen decomposition in compute_spd_sqrt_inv. Max allowed: %d. Consider using diagonal constraints or a method suitable for large sparse matrices.", 
                   A.n_rows, A.n_cols, max_dense_dim);
    }
    // convert to dense; assumed reasonably small if we reach here
    mat M(A);
    vec eigval;
    mat eigvec;
    bool success = eig_sym(eigval, eigvec, M);
    if (!success) {
        Rcpp::stop("Dense eigen decomposition failed in compute_spd_sqrt_inv.");
    }

    uvec keep = find(eigval > tol);
    if (keep.n_elem == 0) {
         Rcpp::warning("No eigenvalues > tol found in compute_spd_sqrt_inv.");
         Ahalf = sp_mat(A.n_rows, A.n_cols);
         Ahalf_inv = sp_mat(A.n_rows, A.n_cols);
         return;
    }
    eigval = eigval(keep);
    eigvec = eigvec.cols(keep);

    vec sqrt_eigval = sqrt(eigval);
    vec inv_sqrt_eigval = 1.0 / sqrt_eigval;

    mat half = eigvec * diagmat(sqrt_eigval) * eigvec.t();
    mat half_inv = eigvec * diagmat(inv_sqrt_eigval) * eigvec.t();

    // Convert back to sparse, maybe check tolerance?
    Ahalf = sp_mat(half);
    Ahalf_inv = sp_mat(half_inv);
}

// operator for Y = R½ Xᵀ Q X R½  (size p×p)
// Note: Assumes X is dense arma::mat
class RightOpProd
{
    const arma::mat&     m_X; // Use reference to avoid copy
    const arma::sp_mat&  m_Q;
    const arma::sp_mat&  m_Rhalf;
    const arma::sp_mat   m_RhalfT; // Keep transpose internally
public:
    RightOpProd(const arma::mat& X_, const arma::sp_mat& Q_, const arma::sp_mat& Rhalf_)
        : m_X(X_), m_Q(Q_), m_Rhalf(Rhalf_), m_RhalfT(Rhalf_.t()) {}

    int rows() const { return m_Rhalf.n_cols; } // Dimension of the operator (p)
    int cols() const { return m_Rhalf.n_cols; } // Dimension of the operator (p)

    // y = Op * x
    void perform_op(const double* x_in, double* y_out) const
    {
        // Map input array to Armadillo vector without copying
        arma::vec v(const_cast<double*>(x_in), cols(), false, true);

        // Perform the matrix-vector products step-by-step
        // No intermediate matrices formed
        arma::vec tmp = m_Rhalf * v;        // R(1/2) * v --> p-dim vector
        tmp = m_X * tmp;                // X * (R(1/2)v) --> n-dim vector
        tmp = m_Q * tmp;                // Q * (X R(1/2) v) --> n-dim vector
        tmp = m_X.t() * tmp;            // X' * (Q X R(1/2) v) --> p-dim vector
        tmp = m_RhalfT * tmp;           // R(1/2)' * (X' Q X R(1/2) v) --> p-dim vector

        // Copy the result into the output array
        std::memcpy(y_out, tmp.memptr(), sizeof(double) * size_t(cols()));
    }
};

// operator for Y = Q½ X R Xᵀ Q½  (size n×n)
// Note: Assumes X is dense arma::mat
class LeftOpProd
{
    const arma::mat&     m_X;
    const arma::sp_mat&  m_Qhalf;
    const arma::sp_mat   m_QhalfT;
    const arma::sp_mat&  m_R;
public:
    LeftOpProd(const arma::mat& X_, const arma::sp_mat& Qhalf_, const arma::sp_mat& R_)
        : m_X(X_), m_Qhalf(Qhalf_), m_QhalfT(Qhalf_.t()), m_R(R_) {}

    int rows() const { return m_Qhalf.n_cols; } // Dimension of the operator (n)
    int cols() const { return m_Qhalf.n_cols; } // Dimension of the operator (n)

    // y = Op * x
    void perform_op(const double* x_in, double* y_out) const
    {
        arma::vec u(const_cast<double*>(x_in), cols(), false, true);

        arma::vec tmp = m_Qhalf * u;        // Q(1/2) * u --> n-dim
        tmp = m_X.t() * tmp;            // X' * (Q(1/2)u) --> p-dim
        tmp = m_R * tmp;                // R * (X' Q(1/2) u) --> p-dim
        tmp = m_X * tmp;                // X * (R X' Q(1/2) u) --> n-dim
        tmp = m_QhalfT * tmp;           // Q(1/2)' * (X R X' Q(1/2) u) --> n-dim

        std::memcpy(y_out, tmp.memptr(), sizeof(double) * size_t(cols()));
    }
};

// ──────────────────────────────────────────────────────────────
// Main routine using Spectra
// ──────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List gmd_fast_cpp(const arma::mat& X,
                  const arma::sp_mat& Q,
                  const arma::sp_mat& R,
                  int k,
                  double tol = 1e-9,
                  int maxit = 1000,
                  int seed = 1234)
{
    //RNGScope __rngScope; // Manage R's RNG state if using R::rnorm etc.
    arma::arma_rng::set_seed(seed); // Set Armadillo's RNG seed

    const int n = X.n_rows;
    const int p = X.n_cols;
    const double eigen_tol = 1e-12; // Tolerance for filtering eigenvalues

    // ---- 1. Compute/Validate Q½, Q⁻½, R½, R⁻½ ----
    sp_mat Qhalf, Qhalf_inv, Rhalf, Rhalf_inv;

    if (Q.is_diagmat()) {
        compute_diag_sqrt_inv(Q, Qhalf, Qhalf_inv, eigen_tol);
    } else {
        compute_spd_sqrt_inv(Q, Qhalf, Qhalf_inv, eigen_tol);
    }

    if (R.is_diagmat()) {
        compute_diag_sqrt_inv(R, Rhalf, Rhalf_inv, eigen_tol);
    } else {
        compute_spd_sqrt_inv(R, Rhalf, Rhalf_inv, eigen_tol);
    }

    // ---- 2. Decide side and setup Spectra solver ----
    bool right_side = (p <= n); // Solve p×p unless p > n
    int dim_op = right_side ? p : n;
    int ncv = std::min(dim_op, std::max(2 * k + 1, 20)); // Default ncv rule for Spectra
    if (ncv <= k) {
      Rcpp::warning("Adjusting ncv (%d) to be > k (%d) for Spectra solver.", ncv, k);
      ncv = std::min(dim_op, k + 1); // Minimal valid ncv
    }
    if (ncv > dim_op) {
        ncv = dim_op; // Cannot exceed dimension
    }

    arma::vec eigval;
    arma::mat W; // Eigenvectors (columns)
    int nconv = 0;

    Rcpp::Rcout << "Starting Spectra eigen solver (k=" << k << ", ncv=" << ncv << ")..." << std::endl;

    if (right_side) {
        // Define the matrix-vector operation for the right-side problem
        RightOpProd op(X, Q, Rhalf);

        // Construct the Spectra solver - pass pointer to op
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, RightOpProd> eigs(&op, k, ncv);
        
        eigs.init(); // Initialize
        try {
           // Compute eigenvalues - SortRule is template param, not arg here
           nconv = eigs.compute(maxit, tol); 
        } catch (const std::exception& e) {
            Rcpp::stop("Spectra (right side) computation failed: %s", e.what());
        } 

        if (eigs.info() != Spectra::SUCCESSFUL) {
            Rcpp::warning("Spectra eigen solver did not converge successfully (right side). Info: %d", (int)eigs.info());
            // Continue with potentially fewer or less accurate eigenvalues?
            // For now, let's proceed but user should be warned.
        }
        
        // Convert Eigen results to Armadillo types via Rcpp::wrap then Rcpp::as
        eigval = Rcpp::as<arma::vec>(Rcpp::wrap(eigs.eigenvalues()));
        W = Rcpp::as<arma::mat>(Rcpp::wrap(eigs.eigenvectors()));
        
    } else { // Left side (n < p)
        // Define the matrix-vector operation for the left-side problem
        LeftOpProd op(X, Qhalf, R);

        // Construct the Spectra solver - pass pointer to op
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, LeftOpProd> eigs(&op, k, ncv);
        
        eigs.init();
        try {
            // Compute eigenvalues - SortRule is template param, not arg here
            nconv = eigs.compute(maxit, tol);
        } catch (const std::exception& e) {
            Rcpp::stop("Spectra (left side) computation failed: %s", e.what());
        }

        if (eigs.info() != Spectra::SUCCESSFUL) {
            Rcpp::warning("Spectra eigen solver did not converge successfully (left side). Info: %d", (int)eigs.info());
        }

        // Convert Eigen results to Armadillo types via Rcpp::wrap then Rcpp::as
        eigval = Rcpp::as<arma::vec>(Rcpp::wrap(eigs.eigenvalues()));
        W = Rcpp::as<arma::mat>(Rcpp::wrap(eigs.eigenvectors()));
    }

    Rcpp::Rcout << "Spectra finished. Converged components: " << nconv << std::endl;

    // ---- 3. Filter eigenvalues and back-transform to U, V, d ----
    uvec keep_idx = find(eigval > eigen_tol);
    int k_found = keep_idx.n_elem;
    
    if (k_found == 0) {
        Rcpp::warning("No positive eigenvalues found > tolerance (%g). Returning empty solution.", eigen_tol);
        return List::create(_["d"] = arma::vec(),
                            _["u"] = arma::mat(n, 0),
                            _["v"] = arma::mat(p, 0),
                            _["k"] = 0);
    }
    
    if (k_found < k) {
         Rcpp::warning("Found only %d eigenvalues > tolerance (%g), less than requested k=%d.", k_found, eigen_tol, k);
    }
    
    // Keep only valid eigenvalues and corresponding eigenvectors
    eigval = eigval(keep_idx);
    W = W.cols(keep_idx);

    arma::vec d = sqrt(eigval); // Generalised singular values
    arma::mat U, V;
    
    if (right_side) {
        // V = R^(-1/2) * W (eigenvectors of right-side problem)
        V = Rhalf_inv * W; // V should be R-orthonormal (V' R V = I)
        
        // U = X V D^(-1), then Q-normalize columns
        // More stable: U_unnorm = X * V; Normalize U_i = U_unnorm[,i] / sqrt(U_i' Q U_i)
        arma::mat U_unnorm = X * V; // n x k_found
        U.set_size(n, k_found);
        for (int i = 0; i < k_found; ++i) {
            arma::vec u_i = U_unnorm.col(i);
            double norm_factor_sq = arma::as_scalar(u_i.t() * Q * u_i);
            if (norm_factor_sq > eigen_tol) {
                U.col(i) = u_i / sqrt(norm_factor_sq);
            } else {
                // Handle near-zero norm? Set to zero or keep as is?
                // Setting to zero might be safer if norm is truly tiny.
                U.col(i).zeros(); 
                Rcpp::warning("Near-zero norm encountered when Q-normalizing U vector %d.", i + 1);
            }
        }
    } else { // Left side
        // U = Q^(-1/2) * W (eigenvectors of left-side problem)
        U = Qhalf_inv * W; // U should be Q-orthonormal (U' Q U = I)
        
        // V = X' U D^(-1), then A-normalize (R-normalize) columns
        // More stable: V_unnorm = X.t() * U; Normalize V_i = V_unnorm[,i] / sqrt(V_i' R V_i)
        arma::mat V_unnorm = X.t() * U; // p x k_found
        V.set_size(p, k_found);
        for (int i = 0; i < k_found; ++i) {
            arma::vec v_i = V_unnorm.col(i);
            double norm_factor_sq = arma::as_scalar(v_i.t() * R * v_i);
             if (norm_factor_sq > eigen_tol) {
                V.col(i) = v_i / sqrt(norm_factor_sq);
            } else {
                V.col(i).zeros();
                Rcpp::warning("Near-zero norm encountered when R-normalizing V vector %d.", i + 1);
            }
        }
    }

    return List::create(_["d"] = d,  // Singular values (k_found)
                        _["u"] = U,  // ou (n x k_found)
                        _["v"] = V,  // ov (p x k_found)
                        _["k"] = k_found // Actual number of components found > tol
                       );
} 