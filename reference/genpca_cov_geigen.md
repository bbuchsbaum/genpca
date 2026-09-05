# Generalized eigenvalue-based covariance GPCA (internal)

Solves the generalized eigenproblem projected onto the retained range of
R. This is the original implementation that was in gpca.R.

## Usage

``` r
genpca_cov_geigen(
  C,
  R = NULL,
  ncomp = NULL,
  constraints_remedy = c("error", "ridge", "clip", "identity"),
  rank_rtol = 1e-06,
  metric_rtol = .metric_rtol_default(),
  verbose = FALSE
)
```
