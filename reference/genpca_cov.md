# Generalized PCA on a covariance matrix (GMD form)

Performs Generalized PCA directly on a pre-computed covariance matrix
`C = X'MX` with a single variable-side metric `R`, following Allen et
al.'s GMD: the eigendecomposition of \\R^{1/2} C R^{1/2}\\ mapped back
with \\V = R^{-1/2} Z\\, so that \\V'RV = I\\. With `C = X'MX` and
`R = A` this matches
[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)`(X, M = M, A = A)`
exactly. This is useful when you already have `C` or when `X` is too
large to store but `C` is manageable.

## Usage

``` r
genpca_cov(
  C,
  R = NULL,
  ncomp = NULL,
  method = c("gmd", "geigen"),
  constraints_remedy = c("error", "ridge", "clip", "identity"),
  rank_rtol = 1e-06,
  metric_rtol = .metric_rtol_default(),
  tol = NULL,
  verbose = FALSE
)
```

## Arguments

- C:

  A p x p symmetric positive semi-definite covariance matrix, typically
  `C = X'MX`. Asymmetry beyond roundoff and indefiniteness beyond
  `metric_rtol` are errors.

- R:

  Variable-side constraint/metric. Can be:

  - NULL: identity matrix (standard PCA on C)

  - a numeric vector of length p: diagonal weights (must be
    non-negative)

  - a p x p symmetric PSD matrix: general metric/smoothing/structure
    penalties

- ncomp:

  Number of components to return. Default is all positive eigenvalues.

- method:

  Deprecated. `"gmd"` (default) is this function; `"geigen"` forwards to
  [`geigen_cov`](https://bbuchsbaum.github.io/genpca/reference/geigen_cov.md)
  with a warning.

- constraints_remedy:

  Deprecated here (GMD requires PSD input and stops otherwise);
  forwarded to
  [`geigen_cov`](https://bbuchsbaum.github.io/genpca/reference/geigen_cov.md)
  when `method = "geigen"`.

- rank_rtol:

  Relative cutoff for component acceptance on the singular-value scale
  (components with `d_j <= rank_rtol * d_1` are dropped). Default 1e-6.

- metric_rtol:

  Relative tolerance for validating `C` and `R` and for detecting the
  numerical null space in an eigendecomposition of a general `R`. Every
  strictly positive diagonal weight is retained without a rank
  approximation. Default `sqrt(.Machine$double.eps)`.

- tol:

  Deprecated; use `rank_rtol` and `metric_rtol`.

- verbose:

  Logical. If TRUE, print progress messages. Default FALSE.

## Value

A plain list (**not** a multivarious `bi_projector`) with components:

- v:

  p x k matrix of loadings (R-orthonormal eigenvectors)

- d:

  Singular values (square root of eigenvalues lambda)

- lambda:

  Eigenvalues (variances under the R-metric)

- k:

  Number of components returned

- propv:

  Proportion of variance explained by each component (total variance is
  \\\mathrm{tr}(CR)\\, Allen et al. Corollary 5)

- cumv:

  Cumulative proportion of variance explained

- R_rank:

  Rank of the constraint matrix R

- method:

  `"gmd"`

Because this is a plain list rather than a `bi_projector`, the
`multivarious` generics `scores()`,
[`components()`](https://bbuchsbaum.github.io/multivarious/reference/components.html),
and
[`reconstruct()`](https://bbuchsbaum.github.io/multivarious/reference/reconstruct.html)
do not apply to it; index `$v`/`$d` directly, or use
[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) when
you need the full projector interface on a data matrix rather than a
pre-computed covariance matrix.

## Details

The generalized eigenproblem \\C v = \lambda R v\\ is a different
estimator (it maximises \\v'Cv\\ subject to \\v'Rv = 1\\, which is
generally gives different components from the GMD) and lives in its own
function,
[`geigen_cov`](https://bbuchsbaum.github.io/genpca/reference/geigen_cov.md).
`method = "geigen"` is accepted here for one release and forwards to it
with a deprecation warning.

## References

Allen, G. I., Grosenick, L., & Taylor, J. (2014). A Generalized
Least-Squares Matrix Decomposition. Journal of the American Statistical
Association, 109(505), 145-159.

## See also

[`geigen_cov`](https://bbuchsbaum.github.io/genpca/reference/geigen_cov.md)
for the generalized eigenproblem \\C v = \lambda R v\\,
[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) for
the two-sided GPCA on data matrices,
[`genpls`](https://bbuchsbaum.github.io/genpca/reference/genpls.md) for
generalized partial least squares

## Examples

``` r
# Standard PCA on a covariance (no constraint)
C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
fit0 <- genpca_cov(C, R=NULL, ncomp=3)
print(fit0$d[1:3])       # first 3 singular values
#> [1] 2.0562689 0.4926162 0.2796596
print(fit0$propv[1:3])   # variance explained by first 3 components
#> [1] 0.92461872 0.05306648 0.01710261

# Equivalence with genpca()
set.seed(123)
X <- matrix(rnorm(50 * 10), 50, 10)
M_diag <- runif(50, 0.5, 1.5)  # row weights
A_diag <- runif(10, 0.5, 2)    # column weights
fit_gpca <- genpca(X, M = M_diag, A = A_diag, ncomp = 5,
                   preproc = multivarious::pass())
C <- crossprod(X, diag(M_diag) %*% X)  # C = X'MX
fit_cov <- genpca_cov(C, R = A_diag, ncomp = 5)
all.equal(fit_gpca$sdev, fit_cov$d, tolerance = 1e-10)
#> [1] TRUE

# Variable weights via a diagonal metric (iris covariance, 4 variables)
C_iris <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
w <- c(1, 1, 0.5, 2)
fitW <- genpca_cov(C_iris, R = w, ncomp=3)
print(fitW$d[1:3])
#> [1] 1.7972881 0.4903437 0.3173048
```
