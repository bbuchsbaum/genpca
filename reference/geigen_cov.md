# Generalized eigenproblem on a covariance matrix

Maximises \\v'Cv\\ subject to \\v'Rv = 1\\ and the additional constraint
that \\v\\ lies in the retained range of `R`. Successive components are
\\R\\-orthogonal. If \\P\\ is the orthogonal projector onto that range,
the returned vectors satisfy \\P C v = \lambda R v\\ and \\V'RV = I\\.
For a full-rank `R` this is the usual equation \\C v = \lambda R v\\. It
also holds for a singular `R` when `C` maps its retained range into
itself. Otherwise the component of \\C v\\ outside the retained range
need not vanish.

## Usage

``` r
geigen_cov(
  C,
  R = NULL,
  ncomp = NULL,
  constraints_remedy = c("error", "ridge", "clip", "identity"),
  rank_rtol = 1e-06,
  metric_rtol = .metric_rtol_default(),
  verbose = FALSE
)
```

## Arguments

- C:

  A p x p symmetric matrix. Asymmetry beyond roundoff is an error; an
  indefinite `C` is allowed (the problem is still defined) and only
  produces a warning.

- R:

  Variable-side constraint/metric. Can be:

  - NULL: identity matrix (standard PCA on C)

  - a numeric vector of length p: diagonal weights (must be
    non-negative)

  - a p x p symmetric PSD matrix: general metric/smoothing/structure
    penalties

- ncomp:

  Number of components to return. Default is all positive eigenvalues.

- constraints_remedy:

  What to do with an indefinite `R`: `"error"` (default), `"ridge"`,
  `"clip"` or `"identity"`; a repair emits a `genpca_metric_repaired`
  warning. See
  [`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md).

- rank_rtol:

  Relative cutoff for component acceptance on the singular-value scale
  (components with `d_j <= rank_rtol * d_1` are dropped). Default 1e-6.

- metric_rtol:

  Relative tolerance for validating `C` and `R` and for detecting the
  numerical null space in an eigendecomposition of a general `R`. Every
  strictly positive diagonal weight is retained without a rank
  approximation. Default `sqrt(.Machine$double.eps)`.

- verbose:

  Logical. If TRUE, print progress messages. Default FALSE.

## Value

A plain list with the same components as
[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
(`v`, `d`, `lambda`, `k`, `propv`, `cumv`, `R_rank`) and
`method = "geigen"`. `propv` is relative to \\\mathrm{tr}(R^{-1/2} C
R^{-1/2})\\ on the range of `R`.

## Details

This is a different estimator from the GMD of
[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md),
which uses \\R^{1/2} C R^{1/2}\\. If `C` and `R` commute, their common
eigenvectors can be ordered differently: the GMD weights variances by
metric eigenvalues, whereas this estimator divides by them. With
`R = c * I`, the directions and their ordering agree, but the eigenvalue
scales differ unless `c = 1`.

`C` is validated for symmetry but may be indefinite (the generalized
eigenproblem is still defined); a warning is issued when its minimum
eigenvalue is below `-metric_rtol * scale`. `R` must be positive
semi-definite; an indefinite `R` is subject to `constraints_remedy`.

## See also

[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)

## Examples

``` r
C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
w <- c(1, 1, 0.5, 2)
fit_gmd <- genpca_cov(C, R = w, ncomp = 2)
fit_geigen <- geigen_cov(C, R = w, ncomp = 2)
# different estimators: the singular values generally differ
rbind(gmd = fit_gmd$d, geigen = fit_geigen$d)
#>            [,1]      [,2]
#> gmd    1.797288 0.4903437
#> geigen 2.658375 0.4966284

# With singular R, the equation is projected onto its retained range
C <- matrix(c(2, 1, 1, 2), 2)
R <- diag(c(1, 0))
fit <- geigen_cov(C, R, ncomp = 1)
P <- diag(c(1, 0))
P %*% C %*% fit$v - (R %*% fit$v) * fit$lambda
#> 2 x 1 Matrix of class "dgeMatrix"
#>      [,1]
#> [1,]    0
#> [2,]    0
```
