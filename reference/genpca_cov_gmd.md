# GMD-based covariance GPCA (internal)

Implements Allen et al.'s GMD approach for covariance matrices. Computes
eigendecomposition of \\R^{1/2} C R^{1/2}\\ and maps back.

## Usage

``` r
genpca_cov_gmd(
  C,
  R = NULL,
  ncomp = NULL,
  rank_rtol = 1e-06,
  metric_rtol = .metric_rtol_default(),
  verbose = FALSE
)
```
