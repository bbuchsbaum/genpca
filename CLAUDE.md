# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

`genpca` is an R package implementing Generalized Principal Component
Analysis and related matrix decompositions with row metric M and column
metric A, following Allen, Grosenick & Taylor (2014). When M = I and A =
I, methods reduce to standard PCA/PLS.

## Development Commands

``` bash
# Load for interactive development
R -e "devtools::load_all()"

# Generate documentation from roxygen comments
R -e "devtools::document()"

# Run all tests
R -e "testthat::test_package('genpca')"

# Run single test file
R -e "testthat::test_file('tests/testthat/test_gpca.R')"

# Full package check
R -e "devtools::check()"

# Build and check (CRAN-style)
R CMD build .
R CMD check genpca_*.tar.gz
```

## Architecture

### Core Statistical Methods

| File             | Function                                                                                                                                       | Purpose                                       |
|------------------|------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------|
| `gpca.R`         | [`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)                                                                          | Generalized PCA with row/column metrics       |
| `gpca_cov.R`     | [`genpca_cov()`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)                                                                  | Covariance-based GPCA (pre-computed C = X’MX) |
| `genpls.R`       | [`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md), [`genplsc()`](https://bbuchsbaum.github.io/genpca/reference/genplsc.md) | Two-block generalized PLS-SVD                 |
| `gplssvd_op.R`   | [`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md)                                                                  | Memory-efficient PLS via implicit operators   |
| `gep_subspace.R` |                                                                                                                                                | Generalized eigenvalue problem methods        |
| `sfpca.R`        |                                                                                                                                                | Sparse functional PCA                         |
| `rpls.R`         |                                                                                                                                                | Regularized partial least squares             |

### Computational Backends

[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
supports multiple methods via the `method` parameter: - `"eigen"`:
Direct eigendecomposition (small-to-medium) - `"spectra"`: Matrix-free
C++ via RSpectra (large/sparse) - `"deflation"`: Sequential extraction
(memory-constrained)

### Key Internal Components

**Constraint Handling** (`gpca.R` + `constraints_utils.R`): -
`prep_constraints()`: Validates and prepares A/M matrices -
[`ensure_spd()`](https://bbuchsbaum.github.io/genpca/reference/ensure_spd.md):
Enforces positive semi-definiteness via Gershgorin shift or nearPD -
[`is_spd()`](https://bbuchsbaum.github.io/genpca/reference/is_spd.md):
Tests SPD via Cholesky decomposition - Remediation strategies:
`"error"`, `"ridge"`, `"clip"`, `"identity"`

**Metric Operators** (`weight_operators.R` + `internal_ops.R`): -
`.metric_operators()`: Creates closures for W, W^{1/2}, W^{-1/2}
multiplication - Optimized paths for diagonal matrices (element-wise
ops) - Cholesky-based paths for general SPD matrices - Used by both GPCA
and GPLSSVD to avoid materializing whitened matrices

**PLS Operators** (`gplssvd_op.R`): - `.build_pls_operator()`:
Constructs S_mv and ST_mv functions for implicit SVD - Computes SVD of S
= t(Xe) %\*% Ye without forming Xe = Mx^{1/2} X Ax^{1/2}

### C++ Layer (src/)

- `gpca.cpp`: Core GPCA algorithms
- `gmd_fast.cpp`: Fast generalized matrix decomposition
- Uses RcppArmadillo and RcppEigen for linear algebra

### Class Hierarchy

Returns S3 objects extending `multivarious` framework: - `genpca` →
`bi_projector` → `projector` - `genpls` → `cross_projector` →
`projector`

Key methods: `scores()`,
[`components()`](https://bbuchsbaum.github.io/multivarious/reference/components.html),
`project()`,
[`reconstruct()`](https://bbuchsbaum.github.io/multivarious/reference/reconstruct.html),
[`transfer()`](https://bbuchsbaum.github.io/multivarious/reference/transfer.html)

## Testing Conventions

Tests verify algorithm equivalence across backends: - Eigen vs Spectra:
sdev within 1e-6, scores within 1e-5 (up to sign) - Deflation vs Eigen:
sdev within 1e-4, subspace agreement via principal angles - CCA vs
genpls equivalence tests - genpls vs gplssvd_op consistency tests

## Key Dependencies

- `Matrix`: Sparse matrix representations (critical for efficiency)
- `multivarious`: Framework for projector classes and preprocessing
- `RSpectra`: Iterative eigensolvers for large problems
- `Rcpp`/`RcppArmadillo`/`RcppEigen`: C++ backend
