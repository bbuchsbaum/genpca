# genpca Package Test Failure Bug Reports

## Report 1: NA Eigenvalue Errors in prep_constraints Function

### Summary
The genpca package is experiencing failures in the `prep_constraints` function where eigenvalue computations are either failing completely (returning NA) or returning negative eigenvalues for matrices that should be positive semi-definite (PSD). This affects multiple test cases involving sparse constraint matrices.

### Root Cause Analysis

#### Primary Issues Identified:

1. **RSpectra eigenvalue computation failures**: The `prep_constraints` function uses `RSpectra::eigs_sym(A, k=1, which="SA")` to find the smallest eigenvalue. This is failing in several scenarios, particularly with sparse matrices created by the `neighborweights` package functions.

2. **Non-PSD matrices from external dependencies**: Test cases are creating constraint matrices using `neighborweights` package functions that are not guaranteed to be positive semi-definite:
   - `neighborweights:::adjacency.neighbor_graph()` creates adjacency matrices 
   - `temporal_adjacency()` followed by `t(A) %*% A` operations
   - These matrices may have negative eigenvalues or numerical issues

3. **Missing dependency**: The `neighborweights` package is used extensively in tests but is not listed in DESCRIPTION as a dependency (Suggests or Imports).

4. **Eigenvalue tolerance issues**: The current tolerance of `1e-8` may be too strict for numerically challenging sparse matrices.

### Specific Failing Tests Analysis:

1. **test_gpca.R:82** - "gen_pca with sparse column and row constraints works"
   - M matrix created from `neighborweights` adjacency graph with `diag(M) <- 1.5`
   - Error: smallest eigenvalue = -0.0283834767954287

2. **test_gpca.R:129** - "can run genpca with sparse weighting matrix" 
   - Large temporal adjacency matrix (10000x10000) 
   - Error: NA (eigenvalue computation failed)

3. **test_gpca.R:144** - "can run genpca on a largeish matrix with deflation"
   - Temporal adjacency matrix `t(A) %*% A` for 500x500 matrix
   - Error: NA (eigenvalue computation failed)

4. **test_gpca.R:243** - "genpca with spatial adjacency recovers a smooth temporal blob"
   - Sparse adjacency matrix from spatial grid
   - Error: smallest eigenvalue = -3.75877048314363

### Impact Assessment

**Severity**: High - Core functionality broken for realistic use cases
**Scope**: 
- Affects users using spatial/temporal constraints
- Breaks large-scale applications with sparse constraint matrices
- 4 critical test failures in main functionality

### Proposed Solutions

1. **Improve eigenvalue computation robustness**
   - Add better fallback strategies when RSpectra fails
   - Implement more robust matrix conditioning checks
   - Use multiple eigenvalue computation methods

2. **Add matrix remediation options**
   - Implement regularization for near-singular matrices
   - Add automatic matrix conditioning (e.g., ridge penalty)
   - Provide better error messages with suggested fixes

3. **Fix dependency issues**
   - Add `neighborweights` to Suggests in DESCRIPTION
   - Or replace with internal functions for test matrix generation

---

## Report 2: dsyMatrix C++ Compatibility Error in gmd_fast_cpp

### Summary
The `gmd_fast_cpp` function in the genpca package fails when passed dense symmetric matrices (dsyMatrix objects) from R's Matrix package, throwing the error "dsyMatrix is not supported". This error occurs during the automatic type conversion from R to C++ through RcppArmadillo.

### Technical Root Cause

The issue stems from incompatible matrix type handling between R's Matrix package and RcppArmadillo's automatic type conversion system:

1. **Matrix Type Creation**: When `Matrix::Matrix()` is called on a dense, symmetric, positive-definite matrix, it automatically creates a `dsyMatrix` object - a specialized class for dense symmetric matrices.

2. **C++ Function Signature**: The `gmd_fast_cpp` function expects `arma::sp_mat` (sparse matrix) types for the Q and R constraint matrices:
   ```cpp
   List gmd_fast_cpp(const arma::mat& X,
                     const arma::sp_mat& Q,    // <-- expects sparse matrix
                     const arma::sp_mat& R,    // <-- expects sparse matrix
                     int k, ...)
   ```

3. **Conversion Failure**: RcppArmadillo's automatic conversion system cannot convert a `dsyMatrix` object to `arma::sp_mat`, resulting in the "dsyMatrix is not supported" error.

### Failing Tests

1. **test_gpca.R:309** - "gmd_fast_cpp matches genpca (use_cpp=TRUE) for p <= n, dense constraints"
2. **test_gpca.R:329** - "gmd_fast_cpp matches genpca (use_cpp=TRUE) for p > n, dense constraints"

### Impact Assessment

**Severity**: High
- Core C++ functionality is completely broken for dense constraint matrices
- Multiple critical tests are failing
- Users cannot use the fast C++/Spectra implementation with dense constraints

### Proposed Solution

Modify the C++ function to accept `SEXP` parameters and handle type conversion explicitly:

```cpp
// In gmd_fast.cpp
List gmd_fast_cpp(const arma::mat& X,
                  SEXP Q_sexp,          // Accept raw SEXP
                  SEXP R_sexp,          // Accept raw SEXP  
                  int k, double tol = 1e-9, int maxit = 1000, int seed = 1234)
{
    // Convert Q to sparse matrix with proper type handling
    arma::sp_mat Q;
    if (Rf_inherits(Q_sexp, "dsyMatrix") || Rf_inherits(Q_sexp, "dgeMatrix")) {
        // Convert dense Matrix to sparse
        arma::mat Q_dense = Rcpp::as<arma::mat>(Q_sexp);
        Q = arma::sp_mat(Q_dense);
    } else {
        Q = Rcpp::as<arma::sp_mat>(Q_sexp);
    }
    
    // Similar handling for R
    // ... rest of function
}
```

---

## Report 3: Numerical Precision Issues in genpca Package Tests

### Summary
Multiple test failures in the genpca R package are occurring due to numerical precision and algorithmic issues across several core functions. The failures span matrix decomposition algorithms, generalized eigenvalue solvers, sparse functional PCA, and subspace comparison utilities.

### Root Cause Analysis

#### 1. Generalized Eigenvalue Problem Solver (`solve_gep_subspace`)
**Issue:** The Gram matrix orthogonality check is failing with large deviations from expected identity matrices.

**Technical cause:**
- The subspace iteration algorithm is not maintaining proper orthogonality in the S2-inner product
- The orthonormalization step uses standard QR decomposition but should use the generalized inner product
- Regularization parameters may be inadequate for ill-conditioned matrices

**Evidence:**
```
gram matrix expected: diag(3) = [1,0,0; 0,1,0; 0,0,1]
gram matrix actual: [0.375, *, 3.496; *, 0.581, *; *, *, *]
```

#### 2. Sparse Functional PCA (`sfpca`)
**Issue:** Rank-1 matrix recovery is severely inaccurate with large magnitude errors.

**Technical cause:**
- The alternating optimization algorithm is not converging to the true rank-1 decomposition
- Regularization penalties are insufficient for the test case
- Sign ambiguity handling is inadequate

**Evidence:**
```
Expected singular value: 17.4
Actual singular value: 1.3
Error magnitude: 92.5% relative error
```

#### 3. C++ vs R Implementation Mismatches (`gmd_fast_cpp`)
**Issue:** Subspace correlation matrices between C++ and R implementations are not close to permutation matrices.

**Technical cause:**
- Different numerical precision or algorithmic paths between C++ and R implementations
- Possible differences in eigenvalue/singular value solver convergence criteria
- Sign convention inconsistencies between implementations

### Common Patterns Across Failures:

1. **Orthogonality Loss:** Multiple algorithms fail to maintain proper orthogonality constraints
2. **Regularization Inadequacy:** Default regularization parameters appear insufficient for test matrices
3. **Sign Ambiguity:** Inconsistent handling of sign indeterminacy in matrix decompositions
4. **Tolerance Mismatches:** Test tolerances don't match the actual numerical precision of algorithms
5. **Implementation Divergence:** C++ and R implementations produce different results for identical inputs

### Proposed Solutions

1. **Fix Generalized Orthogonalization** (Critical)
   - Implement proper generalized Gram-Schmidt orthogonalization in `solve_gep_subspace`
   - Replace standard QR with S2-orthogonal QR decomposition

2. **Improve Regularization Strategy** (High)
   - Implement adaptive regularization based on matrix condition numbers
   - Add condition number monitoring and warnings
   - Increase default regularization parameters

3. **Enhance SFPCA Algorithm** (High)
   - Review convergence criteria in the alternating optimization
   - Implement better initialization strategies
   - Improve penalty parameter estimation heuristics

4. **Standardize Implementations** (Medium)
   - Ensure C++ and R implementations use identical algorithms and tolerances
   - Add comprehensive cross-validation tests between implementations
   - Standardize sign conventions across all matrix decomposition functions

---

## Report 4: RPLS and Transfer Method Issues

### Summary
The genpca R package is experiencing several critical test failures related to missing imports, method dispatch issues, and dimension mismatches in the rpls (regularized partial least squares) and transfer method implementations.

### Root Cause Analysis

#### 1. Missing Function Imports (Primary Issue)

**Problem**: The package does not import essential functions from `multivarious`:
- `project` function is not imported, causing `Error: could not find function "project"`
- `partial_project` function is not imported, causing `Error: object 'partial_project' not found`

**Evidence**: `/Users/bbuchsbaum/code/genpca/NAMESPACE` only imports limited functions from multivarious

#### 2. Method Dispatch and getS3method Issue

**Problem**: `/Users/bbuchsbaum/code/genpca/R/transfer_wrappers.R` line 43 uses `getS3method()` without importing it from the `utils` package.

#### 3. S3 Method Signature Inconsistency

**Problem**: The `transfer.cross_projector` method has parameter name mismatches with the base `multivarious` implementation.

### Failing Tests

1. **test_rpls.R:176** - "rpls partial projection works if partial cols are requested" - Dimension mismatch
2. **test_gpca.R:29** - "gen_pca with column variances is equivalent to a scaled pca" - Attributes length mismatch
3. Package check warnings about S3 method consistency
4. Note about undefined global function 'getS3method'

### Impact Assessment

**Severity**: HIGH
- Multiple core functionality tests failing
- Package fails basic R CMD check requirements
- Method dispatch broken for transfer operations

### Proposed Solutions

#### 1. Fix NAMESPACE Imports (Critical - Priority 1)

Add missing imports to `/Users/bbuchsbaum/code/genpca/NAMESPACE`:
```r
importFrom(multivarious,project)
importFrom(multivarious,partial_project)
importFrom(utils,getS3method)
```

#### 2. Fix transfer_wrappers.R Method (Critical - Priority 1)

Add to `/Users/bbuchsbaum/code/genpca/R/transfer_wrappers.R`:
```r
#' @importFrom utils getS3method
```

#### 3. Standardize Parameter Names (Medium - Priority 2)

Ensure parameter signature consistency in transfer methods.

---

## Summary Statistics

- **Total Test Failures**: 34
- **Categories**: 4 major issue categories
- **Severity**: All HIGH severity
- **Primary Causes**: 
  - Missing imports and dependencies
  - Matrix type compatibility issues
  - Numerical precision problems
  - Algorithm implementation differences

## Recommended Action Plan

1. **Immediate** (Day 1):
   - Fix NAMESPACE imports
   - Add missing dependencies to DESCRIPTION
   - Fix getS3method import

2. **Short-term** (Week 1):
   - Fix dsyMatrix C++ compatibility
   - Improve eigenvalue computation robustness
   - Fix generalized orthogonalization

3. **Medium-term** (Week 2-3):
   - Standardize C++/R implementations
   - Improve SFPCA convergence
   - Update test tolerances

4. **Long-term** (Month 1):
   - Add comprehensive numerical diagnostics
   - Improve documentation
   - Add user guidance for edge cases