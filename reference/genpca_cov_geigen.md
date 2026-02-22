# Generalized eigenvalue-based covariance GPCA (internal)

Solves the generalized eigenproblem C v = lambda R v directly. This is
the original implementation that was in gpca.R.

## Usage

``` r
genpca_cov_geigen(
  C,
  R = NULL,
  ncomp = NULL,
  constraints_remedy = c("error", "ridge", "clip", "identity"),
  tol = 1e-08,
  verbose = FALSE
)
```
