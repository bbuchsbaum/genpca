# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`genpca` is an R package implementing Generalized Principal Component Analysis and related matrix decompositions. The package provides methods for generalized PCA with row and column constraints, as well as generalized partial least squares techniques.

## Development Commands

### Building and Installing the Package
```bash
# Install package with devtools (as configured in .Rproj)
R -e "devtools::install()"

# Build and check package
R CMD build .
R CMD check genpca_*.tar.gz

# Install dependencies
R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'RcppEigen', 'RSpectra', 'FNN', 'Matrix', 'multivarious', 'assertthat'))"
```

### Running Tests
```bash
# Run all tests
R -e "testthat::test_package('genpca')"

# Run tests for a specific file
R -e "testthat::test_file('tests/testthat/test_gpca.R')"

# Run tests with coverage
R -e "covr::package_coverage()"
```

### Development Workflow
```bash
# Load package for interactive development
R -e "devtools::load_all()"

# Generate documentation from roxygen comments
R -e "devtools::document()"

# Check package
R -e "devtools::check()"
```

## Architecture

### Core Components

1. **Main Statistical Methods** (in R/):
   - `gpca.R`: Generalized PCA implementation with constraint handling
   - `gep_subspace.R`: Generalized eigenvalue problem subspace methods
   - `sfpca.R`: Sparse functional PCA
   - `rpls.R`: Regularized partial least squares
   - `plsutils.R`: Utility functions for PLS methods

2. **Transfer Learning** (in R/):
   - `transfer_methods.R`: Core transfer learning implementations
   - `transfer_wrappers.R`: High-level wrappers for transfer methods

3. **C++ Implementation** (in src/):
   - `gpca.cpp`: Core GPCA algorithms in C++
   - `gmd_fast.cpp`: Fast generalized matrix decomposition
   - Uses RcppArmadillo and RcppEigen for efficient linear algebra

### Key Design Patterns

1. **Constraint Handling**: The `prep_constraints()` function in gpca.R validates and prepares constraint matrices (A and M), ensuring they are positive semi-definite with various remediation strategies.

2. **S3 Class System**: Methods follow R's S3 object system with classes like `genpca`, `cross_projector`, and `projector`.

3. **Integration with multivarious**: The package extends the `multivarious` package's framework for dimensionality reduction methods.

### Testing Strategy

Tests use `testthat` framework and are organized by functionality:
- Basic functionality tests (e.g., `test_gpca.R`)
- Algorithm soundness checks (e.g., `test_genpls_alg_soundness.R`)
- Edge case handling (e.g., `test_deflation_empty.R`)
- Integration tests for transfer learning methods

### Important Notes

- The package uses sparse matrix representations from the Matrix package extensively
- Constraint matrices must be positive semi-definite; the package provides multiple remediation strategies
- The experimental/ directory contains work-in-progress implementations
- Package follows R CMD check standards with proper NAMESPACE management