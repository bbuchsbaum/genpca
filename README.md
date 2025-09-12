# genpca

Generalized PCA and related decompositions for data observed in non‑Euclidean inner‑product spaces. The package supports row/column metrics, sparse constraints, matrix‑free solvers, and generalized PLS (PLS‑SVD) — all integrated with the multivarious ecosystem for preprocessing, projection, and reconstruction.

## Features

- Generalized PCA (`genpca`) with row metric `M` (n × n) and column metric `A` (p × p)
- Multiple backends:
  - `method = "eigen"`: direct eigendecomposition (good for modest sizes)
  - `method = "spectra"`: matrix‑free iterative C++ Spectra solver (scalable)
  - `method = "deflation"`: sequential component extraction (memory‑lean)
- Constraint utilities (symmetric/PSD checks and remedies)
- Generalized PLS / PLS‑SVD:
  - High‑level: `genpls()` (alias `genplsc()`)
  - Low‑level operator: `gplssvd_op()` (no explicit “whitened” matrices)
- Seamless integration with `multivarious` (preprocessing, scores/loadings/components, transfer, reconstruct)

## Installation

```r
# install.packages("devtools")
devtools::install_github("bbuchsbaum/genpca")
```

You’ll also want these runtime dependencies installed:

```r
install.packages(c("Matrix", "RSpectra", "multivarious"))
# Optional for some utilities / tests
install.packages(c("irlba", "knitr", "rmarkdown"))
```

## Quick start

```r
library(genpca)
set.seed(1)
X <- matrix(rnorm(200 * 50), 200, 50)

# Standard PCA via GPCA (identity metrics by default)
fit <- genpca(X, ncomp = 5, preproc = multivarious::center())
fit$sdev                        # singular values
head(multivarious::scores(fit)) # n × k
head(multivarious::components(fit)) # p × k

# Row/column metrics (diagonal weights, sparse structures)
library(Matrix)
M <- Diagonal(nrow(X), x = runif(nrow(X), 0.5, 1.5))
A <- bandSparse(ncol(X), ncol(X), k = -1:1,
                diagonals = list(rep(1.2, ncol(X)), rep(-0.3, ncol(X)-1), rep(-0.3, ncol(X)-1)))
fit_w <- genpca(X, M = M, A = A, ncomp = 5, preproc = multivarious::center())

# Canonical PLS (PLS‑SVD)
Y <- matrix(rnorm(200 * 20), 200, 20)
pls <- genpls(X, Y, ncomp = 3, preproc_x = multivarious::center(), preproc_y = multivarious::center())
pls$d            # canonical correlations (singular values)
dim(pls$vx); dim(pls$vy)  # X/Y weights
```

## Metrics in GPCA

`M` and `A` are symmetric positive semi‑definite (PSD) metrics defining inner products in the observation and variable spaces:

- `||x||_M^2 = x^T M x`,   `d_M(x,y)^2 = (x−y)^T M (x−y)`
- `||v||_A^2 = v^T A v`,   `d_A(v,w)^2 = (v−w)^T A (v−w)`

Common choices: identity (standard PCA), diagonal weights (unequal sampling or feature scaling), covariance/precision‑based metrics (to reflect correlation structure), and sparse kernels/Laplacians (temporal/spatial smoothness). When `M = I` and `A = I`, you recover ordinary PCA.

## Documentation

- Vignettes
  - genpca: Generalized PCA and Related Decompositions (overview)
  - Generalized PLS‑SVD: Explicit Whitening Reference (operator vs explicit whitening)

Build locally:

```r
devtools::build_vignettes()
browseVignettes("genpca")
```

## Testing and guarantees

- Eigen vs Spectra: unit tests assert tight agreement on modest problems (sdev within 1e‑6, scores within 1e‑5 up to sign).
- Deflation vs Eigen: additional tests (n≈60, p≈40, k=8) assert:
  - sdev within 1e‑4,
  - subspace agreement via principal angles,
  - scores/components close after Procrustes alignment (≈ 1e‑3 relative).

Run tests locally:

```r
library(testthat)
library(pkgload)
pkgload::load_all()
testthat::test_dir("tests/testthat")
```

## License

MIT (see LICENSE).

## Contributing

Issues and PRs welcome. Please open a ticket with a minimal example, your R session info, and (if relevant) a pointer to the metric matrices that reproduce the behavior.

